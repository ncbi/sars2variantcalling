import sys
import os

args = sys.argv
toolbox_location = os.path.dirname(os.path.dirname(args[args.index("-s") + 1]))

ref = config["ref"]
snpeff_config=config["snpeff_config"]
vcfEffOnePerLine=config["vcfEffOnePerLine"]
vcf_validator=config["vcf_validator"]

accessions_file_path = 'accs'
products = ["consensus.bam", "consensus.coverage", "consensus.depth", "consensus.fa", "consensus.summary", "ref.bam",
            "ref.depth", "ref.snp_eff.tsv", "ref.snpeff.vcf", "ref.summary", "ref.vcf", "vcfvalidate.done"]


with open(accessions_file_path,'r') as f:
    accessions = [line.strip() for line in f.readlines()]

rule all:
    input: expand("{acc}/{acc}.{product}",acc=accessions,product=products)
    shell: """
ls {input}
"""

rule fastq_dump:
    output: R1="{acc}/{acc}_1.fastq",R2="{acc}/{acc}_2.fastq",single="{acc}/{acc}.fastq"
    log: "LOGS/{acc}.fastq_dump.log"
    threads: 6
    shell: """
if [[ ! -d {wildcards.acc} ]]; then mkdir {wildcards.acc}; fi
fasterq-dump --outdir {wildcards.acc} --split-files --threads {threads} {wildcards.acc} &> {log}
touch {output}
"""

rule trimmed:
    input: R1=rules.fastq_dump.output.R1,R2=rules.fastq_dump.output.R2,single=rules.fastq_dump.output.single,
    output: "{acc}/{acc}_R1.trimmed.fastq","{acc}/{acc}_R1.trimmed.unpaired.fastq","{acc}/{acc}_R2.trimmed.fastq","{acc}/{acc}_R2.trimmed.unpaired.fastq","{acc}/{acc}.trimmed.fastq"
    threads: 6
    log: "LOGS/{acc}.trimmed.log"
    shell: """
if [[ -s {input.R1} || -s {input.R2} ]] ; then 
	java -jar /usr/local/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 -threads {threads} -trimlog {log} \\
		{input.R1} {input.R2} \\
		"{wildcards.acc}/{wildcards.acc}_R1.trimmed.fastq" "{wildcards.acc}/{wildcards.acc}_R1.trimmed.unpaired.fastq" "{wildcards.acc}/{wildcards.acc}_R2.trimmed.fastq" "{wildcards.acc}/{wildcards.acc}_R2.trimmed.unpaired.fastq" \\
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

	java -jar /usr/local/trimmomatic/0.33/trimmomatic-0.33.jar SE -phred33 -threads 6 {input.single} {wildcards.acc}/{wildcards.acc}.trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

if [[ -s input.unpaired ]]; then 
	java -jar /usr/local/trimmomatic/0.33/trimmomatic-0.33.jar SE -phred33 -threads {threads} \\
		-trimlog {log} input.single \\ 
		"{wildcards.acc}/{wildcards.acc}.trimmed.fastq" \\
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

touch {output}
"""

rule bam:
    input: rules.trimmed.output
    output: bam="{acc}/{acc}.ref.bam",bai="{acc}/{acc}.ref.bam.bai",summary="{acc}/{acc}.ref.summary"
    threads: 6
    log: hisat2_log="LOGS/{acc}.hisat2_bam.log",picard_log="LOGS/{acc}.picard_add_groups.log"
    shell: """
tmpbam=$(mktemp {wildcards.acc}.XXX.bam)
( hisat2 --no-spliced-alignment --no-unal -x {ref} -q \\
	-1 "{wildcards.acc}/{wildcards.acc}_R1.trimmed.fastq" -2 "{wildcards.acc}/{wildcards.acc}_R2.trimmed.fastq" \\
	-U "{wildcards.acc}/{wildcards.acc}_R1.trimmed.unpaired.fastq","{wildcards.acc}/{wildcards.acc}_R2.trimmed.unpaired.fastq","{wildcards.acc}/{wildcards.acc}.trimmed.fastq" \\
	--summary-file {output.summary} --threads {threads} | \\
	samtools view -Sb -F256 - | \\
	samtools sort - > $tmpbam) 2>{log.hisat2_log} > {output.bam}

picard AddOrReplaceReadGroups \\
    I=$tmpbam O={output.bam} \\
    RGID=1 RGPl=Illumina RGPU=NA RGSM={wildcards.acc} RGLB=NA >&{log.picard_log}
samtools index {output.bam}
rm $tmpbam
"""

rule call:
    input: bam=rules.bam.output.bam
    output: gvcf="{acc}/{acc}.ref.gvcf"
    log: "LOGS/{acc}.call.log"
    shell: """
gatk HaplotypeCaller -R {ref} -I {input.bam} -O {output.gvcf} --minimum-mapping-quality 10 --ploidy 2 -ERC BP_RESOLUTION >&{log}
"""

rule ref_depth:
    input: rules.call.output.gvcf
    output: "{acc}/{acc}.ref.depth"
    shell: """
vcftools --vcf {input} --extract-FORMAT-info DP --stdout | tail -n +2 > {output}
"""

rule filter_variants:
    input: rules.call.output.gvcf
    output: "{acc}/{acc}.ref.filtered.vcf"
    log: "LOGS/{acc}.filter_variants.log"
    shell: """
tempvcf=$(mktemp {wildcards.acc}.XXX.vcf)
gatk GenotypeGVCFs  -R {ref} -V {input} -O tempvcf

gatk VariantFiltration \
    -R {ref} \
    -V tempvcf \
    -O {output} \
    --filter-name "lowAF" \
    --filter-expression 'vc.getGenotype("{wildcards.acc}").getAD().1.floatValue() / vc.getGenotype("{wildcards.acc}").getDP() < 0.15' \
    --filter-name "lowDP" \
    --filter-expression 'vc.getGenotype("{wildcards.acc}").getDP() < 50' >&{log}
rm $tempvcf
"""

rule norm:
    input: rules.filter_variants.output
    output: "{acc}/{acc}.ref.vcf"
    threads: 1
    log: "LOGS/{acc}.norm.log"
    shell: """
( gatk LeftAlignAndTrimVariants --verbosity ERROR --split-multi-allelics --QUIET \
    -R {ref} \
    -V {input} \
    -O {output} ) >&{log}

grep -vq "^#" {output} || echo "No snps found"
"""

rule genomecov:
    input: bam="{acc}/{acc}.ref.bam"
    output: "{acc}/{acc}.gatk.bam.avg_cov"
    log: "LOGS/{acc}.genomecov.log"
    shell: """
( bedtools genomecov -d -ibam {input} | awk 'BEGIN {{sum=0}}; {{sum+=$3}}; END{{print sum/NR}}' ) 2>{log} > {output}
"""

rule spdi:
    input: rules.norm.output
    output: vcf="{acc}/{acc}.ref.spdi.vcf", summary="{acc}/{acc}.ref.spdi.summary"
    log: "LOGS/{acc}.spdi.log"
    threads: 1
    shell: """
python3 {toolbox_location}/rules/common/SPDI.py --r {ref} --i {input} --o {output.vcf} --s {output.summary}
"""

rule snpeff:
    input: rules.spdi.output.vcf
    output: "{acc}/{acc}.ref.snpeff.vcf"
    log: "LOGS/{acc}.snpeff.log"
    threads: 1
    shell: """
snpeff ann \\
	-nodownload -formatEff -classic -noStats -noLog -quiet -no-upstream -no-downstream \\
	-c {snpeff_config} sars2 -v {input} \\
	2>{log} > {output}
"""

rule tsv:
    input: rules.snpeff.output
    output: "{acc}/{acc}.ref.snp_eff.tsv"
    shell: """
cat {input} | \\
{vcfEffOnePerLine} | \\
snpsift \\
    extractFields - -s "," -e "." \\
    CHROM POS REF ALT \\
    "GEN[0].DP" "GEN[0].AD[1]" \\
    "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" > {output}
"""

rule vcfValidate:
    input: snpvcf = rules.snpeff.output
    output: touch("{acc}/{acc}.vcfvalidate.done")
    log: "LOGS/{acc}.vcf_validate.log"
    threads: 1
    shell:
        """
            {vcf_validator} -i {input.snpvcf} 2>{log}
        """

rule consensus:
    input: rules.norm.output
    output: "{acc}/{acc}.consensus.fa"
    log: "LOGS/{acc}.consensus.log"
    threads: 1
    shell: """
if ! bcftools view {input} -Oz -o {wildcards.acc}/{wildcards.acc}.vcf.gz; then
    echo -n "failed-to-gzip-vcf" 1>&2
    exit 1
fi

if ! bcftools index -f {wildcards.acc}/{wildcards.acc}.vcf.gz; then
    echo -n "failed-to-index-vcf" 1>&2
    exit 1
fi

( bcftools consensus -f {ref} {wildcards.acc}/{wildcards.acc}.vcf.gz | \\
  sed -r "s/^>([[:print:]])*/>{wildcards.acc}_consensus/g" > {output} ) 2>{log}
"""

rule align_consensus:
    input: fastq=rules.trimmed.output,consensus=rules.consensus.output
    output: bam="{acc}/{acc}.consensus.bam",bai="{acc}/{acc}.consensus.bam.bai", summary="{acc}/{acc}.consensus.summary"
    threads: 6
    log: "LOGS/{acc}.ref.bam.log"
    shell: """
hisat2-build {input.consensus} {wildcards.acc}/{wildcards.acc}.index &>/dev/null

( hisat2 --no-spliced-alignment --no-unal -x {ref} -q \\
	-1 "{wildcards.acc}/{wildcards.acc}_R1.trimmed.fastq" -2 "{wildcards.acc}/{wildcards.acc}_R2.trimmed.fastq" \\
	-U "{wildcards.acc}/{wildcards.acc}_R1.trimmed.unpaired.fastq","{wildcards.acc}/{wildcards.acc}_R2.trimmed.unpaired.fastq","{wildcards.acc}/{wildcards.acc}.trimmed.fastq" \\
	--summary-file {output.summary} --threads {threads} | samtools view -Sh -F 256 - | samtools sort - >{output.bam}) 2>{log}
	
samtools index {output.bam}
"""

rule get_coverage:
    input: rules.align_consensus.output.bam
    output: "{acc}/{acc}.consensus.coverage"
    shell: """
samtools coverage {input} > {output}
"""

rule get_depth:
    input: rules.align_consensus.output.bam
    output: "{acc}/{acc}.consensus.depth"
    shell: """
samtools depth -aa -d 0 -m 0 {input} > {output}
"""

