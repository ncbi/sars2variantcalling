import json

products = ["filtered_fq", "initial.consensus.fa", "initial.avg_depth", "final.consensus.fa", "mummer.delta",
            "delta.snps", "delta.snps.vcf", "vcf"]

with open(config['jobs'],'r') as jj:
    tickets = json.load(jj)
    for ticket in tickets:
        parameters = tickets[ticket]
        parameters['model'] = "r941_min_high_g360"
        if 'instrument' in parameters and parameters['instrument'] == 'PromethION':
            parameters['model'] = "r941_prom_high_g360"

accessions = list(tickets.keys())

rule all:
    input: expand("{acc}/{acc}.{product}",acc=accessions,product=products)
    shell: """
ls {input}
"""

rule fastq_dump:
    output: "{acc}/{acc}.fq"
    log: "LOGS/{acc}.fastq-dump.log"
    threads: 6
    shell: """
if [[ ! -d {wildcards.acc} ]]; then mkdir {wildcards.acc}; fi

fastq-dump --stdout --accession {wildcards.acc} > {output} 2> {log}
"""

rule nano_filt:
    input: rules.fastq_dump.output
    output: "{acc}/{acc}.filtered_fq"
    log: "LOGS/{acc}.filter.log"
    threads: 6
    shell: """
cat {input} | NanoFilt -q 10 --headcrop 40 > {output} 2> {log} 
"""


def medaka_consensus_ref(wildcards):
    if {wildcards.stage} == 'initial':
        return config[ref]
    else:
        return f"{wildcards.acc}/{wildcards.acc}.final.consensus.fa"


rule medaka_consensus_initial:
    input: fq=rules.nano_filt.output
    output: consensus="{acc}/{acc}.initial.consensus.fa",gaps="{acc}/{acc}.initial.gaps.bed",
        bam="{acc}/{acc}.initial.ref.bam",bai="{acc}/{acc}.initial.ref.bam.bai",
    params: model=lambda wildcards: tickets[wildcards.acc]['model']
    log: "LOGS/{acc}.initial.medaka-consensus.log"
    threads: 6
    shell: """
medaka_consensus.sh \
    {input.fq} \
    {wildcards.acc} \
    {config[ref]} \
    {params.model} \
    initial \
    NC_045512.2 \
    {config[nanopore_coverage_threshold]} \
    {threads} \
    LOGS
"""

rule medaka_consensus_final:
    input: fq=rules.nano_filt.output,ref=rules.medaka_consensus_initial.output.consensus
    output: consensus="{acc}/{acc}.final.consensus.fa",gaps="{acc}/{acc}.final.gaps.bed",
        bam="{acc}/{acc}.final.ref.bam",bai="{acc}/{acc}.final.ref.bam.bai",
    params: model=lambda wildcards: tickets[wildcards.acc]['model']
    log: "LOGS/{acc}.final.medaka-consensus.log"
    threads: 6
    shell: """
medaka_consensus.sh \
    {input.fq} \
    {wildcards.acc} \
    {input.ref} \
    {params.model} \
    final \
    NC_045512.2 \
    {config[nanopore_coverage_threshold]} \
    {threads} \
    LOGS
"""

rule consensus_analysis:
    input: bam="{acc}/{acc}.initial.ref.bam",bed="{acc}/{acc}.initial.gaps.bed"
    output: avg_depth="{acc}/{acc}.initial.avg_depth",genomecov="{acc}/{acc}.initial.genomecov",
        total_gaps="{acc}/{acc}.initial.total.gaps",mean_stddev_gaps="{acc}/{acc}.initial.mean_stddev_gaps"
    shell: """
depth_avg_stddev=$( bedtools genomecov -ibam {input.bam} -g {config[ref]} -d | awk '{{sum+=$3; sumsq+=$3*$3}} END {{print sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}}' )
echo ${{depth_avg_stddev}} > {output.avg_depth}

total_gaps=$( awk '{{sum += ($3-$2+1) }} END {{print sum}}' {input.bed} ) 
[[ -z "$total_gaps" ]] && total_gaps=0

echo ${{total_gaps}} > {output.total_gaps}
 
bedtools genomecov -ibam {input.bam} -g {config[ref]} -bga > {output.genomecov}

echo -n "${{depth_avg_stddev}} ${{total_gaps}}" > {output.mean_stddev_gaps}
"""

rule nucmer:
    input: consensus=rules.medaka_consensus_final.output.consensus
    output: delta="{acc}/{acc}.mummer.delta"
    params: prefix=lambda wildcards: wildcards.acc
    shell: """
nucmer --prefix={params.prefix} {config[ref]} {input.consensus}
if [[ "{params.prefix}.delta" != "{output.delta}" ]]; then 
    mv {params.prefix}.delta {output.delta} 
fi
"""

rule show_snps:
    input: delta=rules.nucmer.output.delta
    output: snps="{acc}/{acc}.delta.snps"
    shell: """
show-snps -T {input.delta} > {output.snps}
"""

rule mummer2vcf:
    input: snps=rules.show_snps.output.snps
    output: vcf="{acc}/{acc}.delta.snps.vcf"
    shell: """
mummer2vcf.py \
    -s {input.snps} --input-header -t ALL -g {config[ref]} > {output.vcf}.no_header 

cat <<EOF > {output.vcf}
##fileformat=VCFv4.1
##contig=<ID=NC_045512.2,length=29903>
##INFO=<ID=INDEL,Number=1,Type=String,Description="VAR TYPE">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
EOF

cat {output.vcf}.no_header | grep -v ^# |cut -f1-8 >> {output.vcf}
"""

rule mpileup:
    input: bam=rules.medaka_consensus_initial.output.bam
    output: pileup="{acc}/{acc}.initial.ref.bam.pileup"
    log: "LOGS/{acc}.mpileup.log"
    shell: """
samtools mpileup \
    --min-MQ 0 --min-BQ 0 --count-orphans --no-output-ends --no-output-del --no-output-ins --ignore-overlaps -B \
    --fasta-ref {config[ref]} \
    {input.bam} > {output.pileup}
"""

rule nano_annot:
    input: pileup=rules.mpileup.output.pileup,vcf=rules.mummer2vcf.output.vcf,
        mean_stddev_gaps=rules.consensus_analysis.output.mean_stddev_gaps
    output: vcf="{acc}/{acc}.vcf"
    threads: 1
    shell: """
read mean stddev gaps <<< $(cat {input.mean_stddev_gaps}; echo )

bcftools filter --SnpGap 10:INDEL {input.vcf} | grep -v INDEL | vcf-annotate --filter c=2,10 > {input.vcf}.filt_with_snpCluster

set +e
snps=$(grep -c -v "^#" {input.vcf}.filt_with_snpCluster)
if ((snps==0)); then
    echo -n "no-SNPs" 1>&2
    exit 1
fi

set -e
set -o pipefail

(
set -e
set -o pipefail 
grep '^##' {input.vcf}.filt_with_snpCluster
cat << EOF
##INFO=<ID=AVG,Number=1,Type=float,Description="Mean coverage of the run">
##INFO=<ID=STD,Number=1,Type=float,Description="Standard devidation of the coverage">
##INFO=<ID=GAP,Number=1,Type=Integer,Description="Total counts of the gaps from 0 coverage">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw depths">
##INFO=<ID=AD,Number=1,Type=Integer,Description="Depth of ALT">
##INFO=<ID=RUN,Number=1,Type=String,Description="SRA run accession">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

# The main idea is to join (in SQL sense) the snpCluster (VCF) and mpileup on the position column.
# Then output only these columns: 
# from the VCF: CHROM POS REF ALT QUAL INFO
# from the pileup: DEPTH BASES  
join -1 2 -2 2 <(grep -v '^#' {input.vcf}.filt_with_snpCluster | grep -v SnpCluster | sort -k2,2) <(sort -k2,2 {input.pileup}) -o "1.1 1.2 1.4 1.5 1.6 1.7 2.4 2.5" | sort -nk2,2 | while read chrom pos ref alt qual filter depth bases; do
  # some boiler plate to simplify the bases from pileup, we just need counts for ACGTN and indels:
  AD=$(echo "$bases" | sed -r 's/\^]//g' | sed -r 's/[+-][0-9]+/*/g' | tr A-Z a-z | grep -o . | sort | uniq -c | grep -wi $alt | awk '{{print $1}}')

  # samtools mpileup always counts deletions as depth whereas other tools don't.
  # For example samtools depth mimics this with -J option.
  # Subtracting the number of deletions at the position to match the samtools depth/bedtools genomecov behaviour:
  del_count=$(echo "$bases" | tr -dc '*' | wc -c)
  ((depth = depth - del_count))

  info="RUN={wildcards.acc};AVG=$mean;STD=$stddev;GAP=$gaps;DP=$depth;AD=$AD"
  echo $chrom $pos . $ref $alt $qual $filter $info | tr ' ' '\t'
done ) > {output.vcf}.tmp
echo 1 >&2
cat {output.vcf}.tmp | grep -v SnpCluster > {output.vcf}
echo 2 >&2
set +e
snps=$(grep -c -v "^#" {output.vcf})
if ((snps==0)); then
    echo -n "no-SNPs" 1>&2
    exit 1
fi
echo 3 >&2

site3=$( cat {output.vcf}.tmp | awk '{{if ($2 == 28881 || $2 == 28882 || $2 == 28883 ) print }}' | wc | awk '{{print $1}}' )
if [[ $site3 -eq 3 ]]; then
    cat {output.vcf}.tmp | awk '{{if ($2 == 28881 || $2 == 28882 || $2 == 28883 ) print }}' >> {output.vcf}
fi
"""
