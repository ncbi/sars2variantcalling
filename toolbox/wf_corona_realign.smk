ref=config["ref"]
accessions_file_path=config['accs']

with open(accessions_file_path, 'r') as f:
        accessions=[line.strip() for line in f.readlines()]

rule all:
	input: expand("INFO/{acc}.info", acc=accessions), expand("CONSENSUS/{acc}.fasta", acc=accessions), expand("CONSENSUS_BAM/{acc}.bam", acc=accessions), expand("CALL/{acc}.low_cov.bed", acc=accessions), expand("ANN/{acc}.gatk.snpeff.tsv", acc=accessions)

rule info:
	output: "INFO/{acc}.info"
	shell: """
vdb-dump --info {wildcards.acc} > {output}
"""

rule fastq_dump:
	output: "FASTQ/{acc}_1.fastq", "FASTQ/{acc}_2.fastq", "FASTQ/{acc}.fastq"
	log: "LOGS/{acc}.fastq_dump.log"
	threads: 6
	shell: """
fasterq-dump --outdir FASTQ --split-files --threads {threads} {wildcards.acc} &> {log}
touch {output}
"""

rule trimmed:
	input: rules.fastq_dump.output
	output: "TRIMMED/{acc}_R1.trimmed.fastq", "TRIMMED/{acc}_R1.trimmed.unpaired.fastq", "TRIMMED/{acc}_R2.trimmed.fastq", "TRIMMED/{acc}_R2.trimmed.unpaired.fastq", "TRIMMED/{acc}.trimmed.fastq"
	threads: 6
	log: "LOGS/{acc}.trimmed.log"
	shell: """
if [[ -s {input[0]} || -s {input[1]} ]] ; then 
	java -jar /usr/local/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 -threads {threads} -trimlog {log} \\
		{input[0]} {input[1]} \\
		"TRIMMED/{wildcards.acc}_R1.trimmed.fastq" "TRIMMED/{wildcards.acc}_R1.trimmed.unpaired.fastq" "TRIMMED/{wildcards.acc}_R2.trimmed.fastq" "TRIMMED/{wildcards.acc}_R2.trimmed.unpaired.fastq" \\
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

if [[ -s input[3] ]]; then 
	java -jar /usr/local/trimmomatic/0.33/trimmomatic-0.33.jar SE -phred33 -threads {threads} \\
		-trimlog {log} input[3] \\ 
		"TRIMMED/{wildcards.acc}.trimmed.fastq" \\
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

touch {output}
"""

rule bam: 
	input: rules.trimmed.output
	output: bam="BAM/{acc}.ref.bam", bai="BAM/{acc}.ref.bam.bai"
	threads: 6
	log: "LOGS/{acc}.ref.bam.log"
	shell: """
tmpbam=$(mktemp -p BAM)
( /home/zaluninvv/hisat2-2.1.0/hisat2 --no-spliced-alignment --no-unal -x {ref} -q \\
	-1 "TRIMMED/{wildcards.acc}_R1.trimmed.fastq" -2 "TRIMMED/{wildcards.acc}_R2.trimmed.fastq" \\
	-U "TRIMMED/{wildcards.acc}_R1.trimmed.unpaired.fastq","TRIMMED/{wildcards.acc}_R2.trimmed.unpaired.fastq","TRIMMED/{wildcards.acc}.trimmed.fastq" \\
	--summary-file BAM/{wildcards.acc}.hisat2.summary --threads {threads} | samtools view -Sh -F 256 - | samtools sort - > $tmpbam) 2>{log} > {output.bam}
	
/usr/local/picard/2.20.1/bin/picard AddOrReplaceReadGroups I=$tmpbam O={output.bam} RGID=1 RGPl=Illumina RGPU=NA RGSM={wildcards.acc} RGLB=NA
samtools index {output.bam}
rm $tmpbam
"""

rule call:
	input: bam=rules.bam.output.bam, bai=rules.bam.output.bai
	output: vcf="CALL/{acc}.gatk.vcf", bam="CALL/{acc}.gatk.bam"
	shell: """
gatk HaplotypeCaller -R {ref} -I {input.bam} -O {output.vcf} -bamout {output.bam} --minimum-mapping-quality 10 --ploidy 1
"""

rule filter_variants:
	input: rules.call.output.vcf
	output: "FILTERED/{acc}.vcf"
	shell: """
gatk VariantFiltration \
    -R {ref} \
    -V {input} \
    -O {output} \
    --filter-name "lowAF" \
    --filter-expression 'vc.getGenotype("{wildcards.acc}").getAD().1.floatValue() / vc.getGenotype("{wildcards.acc}").getDP() < 0.15' \
    --filter-name "lowDP" \
    --filter-expression 'vc.getGenotype("{wildcards.acc}").getDP() < 50'
"""

rule norm:
	input: rules.filter_variants.output
	output: "NORM/{acc}.vcf"
	threads: 1
	log: "LOGS/{acc}.norm.log"
	shell: """
( gatk LeftAlignAndTrimVariants --verbosity ERROR --split-multi-allelics --QUIET \
    -R {ref} \
    -V {input} \
    -O {output} ) 2>{log}

grep -v -m1 "^#" {output} || echo "No snps found"
"""

rule genomecov:
	input: bam="CALL/{acc}.gatk.bam"
	output: "CALL/{acc}.gatk.bam.avg_cov"
	log: "LOGS/{acc}.genomecov.log"
	shell: """
( bedtools genomecov -d -ibam {input} | awk 'BEGIN {{sum=0}}; {{sum+=$3}}; END{{print sum/NR}}' ) 2>{log} > {output}
"""

rule consensus:
	input: rules.norm.output
	output: "CONSENSUS/{acc}.fasta"
	log: "LOGS/{acc}.consensus.log"
	threads: 1
	shell: """
if ! bcftools view {input} -Oz -o CONSENSUS/{wildcards.acc}.vcf.gz; then
    echo -n "failed-to-gzip-vcf" 1>&2
    exit 1
fi

if ! bcftools index -f CONSENSUS/{wildcards.acc}.vcf.gz; then
    echo -n "failed-to-index-vcf" 1>&2
    exit 1
fi

( bcftools consensus -f {ref} CONSENSUS/{wildcards.acc}.vcf.gz | sed -r "s/^>([[:print:]])*/>{wildcards.acc}_consensus/g" > {output} ) 2>{log}
"""

rule align_consensus:
	input: fastq=rules.trimmed.output, consensus=rules.consensus.output
	output: bam="CONSENSUS_BAM/{acc}.bam", bai="CONSENSUS_BAM/{acc}.bam.bai"
	threads: 6
	log: "LOGS/{acc}.ref.bam.log"
	shell: """
hisat2-build {input.consensus} CONSENSUS_BAM/{wildcards.acc}.index &>/dev/null

( /home/zaluninvv/hisat2-2.1.0/hisat2 --no-spliced-alignment --no-unal -x {ref} -q \\
	-1 "TRIMMED/{wildcards.acc}_R1.trimmed.fastq" -2 "TRIMMED/{wildcards.acc}_R2.trimmed.fastq" \\
	-U "TRIMMED/{wildcards.acc}_R1.trimmed.unpaired.fastq","TRIMMED/{wildcards.acc}_R2.trimmed.unpaired.fastq","TRIMMED/{wildcards.acc}.trimmed.fastq" \\
	--summary-file {wildcards.acc}.hisat2.summary --threads {threads} | samtools view -Sh -F 256 - | samtools sort - >{output.bam}) 2>{log}
	
samtools index {output.bam}
"""

rule get_coverage:
	input: rules.align_consensus.output.bam
	output: "CONSENSUS_BAM/{acc}.consensus.coverage"
	shell: """
samtools coverage {input} > {output}
"""

rule get_depth:
	input: rules.align_consensus.output.bam
	output:  "CONSENSUS_BAM/{acc}.consensus.depth"
	shell: """
samtools depth -aa -d 0 -m 1000000 {input} > {output}
"""

rule snpeff:
	input: rules.norm.output
	output: "ANN/{acc}.vcf"
	log: "LOGS/{acc}.snpeff.log"
	threads: 1
	shell: """
/panfs/traces01.be-md.ncbi.nlm.nih.gov/sra_review/scratch/zaluninvv/Test/snpEff/snpEff/exec/snpeff ann \\
	-nodownload -formatEff -classic -noStats -noLog -quiet -no-upstream -no-downstream -c /home/zaluninvv/snpEff/snpEff.config -dataDir /home/zaluninvv/snpEff/data/ sars2 -v {input} \\
	2>{log} > {output}
"""

rule tsv:
	input: rules.snpeff.output
	output: "ANN/{acc}.gatk.snpeff.tsv"
	shell: """
cat {input} | /usr/local/snpEff/4.2/scripts/vcfEffOnePerLine.pl | /panfs/traces01.be-md.ncbi.nlm.nih.gov/sra_review/scratch/zaluninvv/Test/snpEff/snpEff/exec/snpsift extractFields - -s "," -e "." CHROM POS REF ALT {wildcards.acc} "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" > {output}
"""

rule  low_cov:
	input: rules.call.output.bam
	output: "CALL/{acc}.low_cov.bed"
	shell: """
bedtools genomecov -max 10 -bga -ibam {input} | awk '$4<1' | bedtools merge -i - > {output}
"""
