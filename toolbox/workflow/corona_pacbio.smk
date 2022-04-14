ref = config["ref"]
snpeff_config=config["snpeff_config"]
vcfEffOnePerLine=config["vcfEffOnePerLine"]
vcf_validator=config["vcf_validator"]

accessions_file_path = 'accs'
products = ["fastq", "trimmed.fastq", "bam_stats", "bam_genomecov", "bam_avg_std_depth", "bam_gaps"]

#products = ["fastq", "consensus.bam", "consensus.coverage", "consensus.depth", "consensus.fa", "consensus.summary", "ref.bam",
#            "ref.depth", "ref.snp_eff.tsv", "ref.snpeff.vcf", "ref.summary", "ref.vcf", "vcfvalidate.done"]

with open(accessions_file_path,'r') as f:
    accessions = [line.strip() for line in f.readlines()]

rule all:
    input: expand("{acc}/{acc}.{product}",acc=accessions,product=products)
    shell: """
ls -l {input}
"""

rule fastq_dump:
    output: R="{acc}/{acc}.fastq"
    log: "LOGS/{acc}.fastq_dump.log"
    threads: 6
    shell: """
if [[ ! -d {wildcards.acc} ]]; then mkdir {wildcards.acc}; fi
fasterq-dump --outdir {wildcards.acc} --threads {threads} {wildcards.acc} &> {log}
touch {output}
"""

rule trimmed:
    input: rules.fastq_dump.output.R
    output: "{acc}/{acc}.trimmed.fastq"
    log: "LOGS/{acc}.trimmed.log"
    shell: """
cat {input} | awk '{{if( NR % 4 == 2 || NR % 4 == 0 ) {{print substr ($0, 40, length ($0) - 80) }} else {{print $0}} }}' > {output}
"""

rule bam:
    input: rules.trimmed.output
    output: bam="{acc}/{acc}.bam"
    shell: """
RG="@RG\tID:{wildcards.acc}\tSM:{wildcards.acc}\tLB:{wildcards.acc}\tPL:Pacbio"
minimap2 -ax asm20 --MD -R "$RG" {ref} {input} | samtools sort -@ 2 -o {output.bam}
samtools index {output.bam}
"""

rule bam_stats:
    input: rules.bam.output
    output: "{acc}/{acc}.bam_stats"
    shell: """
samtools stats {input} > {output}
"""

rule genome_cov:
    input: rules.bam.output
    output: avg_std="{acc}/{acc}.bam_avg_std_depth",g_cov="{acc}/{acc}.bam_genomecov"
    shell: """
bedtools genomecov -ibam {input} -g {ref} -d |  awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}}' > {output.avg_std}
bedtools genomecov -ibam {input} -g {ref} -bga > {output.g_cov}
"""

rule gaps:
    input: rules.genome_cov.output.g_cov
    output: "{acc}/{acc}.bam_gaps"
    shell: """
g=$( cat {input} | awk '{{if($4 == 0) sum += ($3-$2) }} END {{print sum}}' )
g_pct=$(( 100 * $g / 29903 ))
echo "gaps = $g, g_pct = $g_pct for run {wildcards.acc} ....."
echo "$g $g_pct" >> {output}
if ((g_pct > 20)); then
    echo -n "not-enough-coverage" 1>&2
    exit 1
fi
"""

