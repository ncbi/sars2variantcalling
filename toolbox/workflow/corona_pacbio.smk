ref = config["ref"]
snpeff_config=config["snpeff_config"]
vcfEffOnePerLine=config["vcfEffOnePerLine"]
vcf_validator=config["vcf_validator"]
snpEff=config["snpEff"]
snpsift=config["snpsift"]
gatk=config["gatk"]

accessions_file_path = 'accs'
products = ["fastq", "trimmed.fastq", "bam_stats", "bam_genomecov",
            "bam_avg_std_depth", "bam_gaps", "missmatch", "ref.depth", "ref.snp_eff.tsv", "ref.snpeff.vcf",
            "vcfvalidate.done"]

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
RG="@RG\\tID:{wildcards.acc}\\tSM:{wildcards.acc}\\tLB:{wildcards.acc}\\tPL:Pacbio"
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
g=$( cat {input} | awk '{{if($4 == 0) {{sum += ($3-$2)}} else {{sum += 0}} }} END {{print sum}}' )
g_pct=$(( 100 * $g / 29903 ))
echo "gaps = $g, g_pct = $g_pct for run {wildcards.acc} ....."
echo "$g $g_pct" >> {output}
if ((g_pct > 20)); then
    touch {wildcards.acc}/low_genome_coverage_gaps_20pct_plus
    echo -n "not-enough-coverage" 1>&2
    exit 1
fi
"""

rule missmatch:
    input: rules.bam_stats.output
    output: "{acc}/{acc}.missmatch"
    shell: """
mismatches=`cat {input} | grep mismatches | grep "NM fields" | awk '{{print $3}}'`
bases_mapped=`cat {input}  | grep cigar | grep "more accurate" | awk '{{print $5}}'`

if ((bases_mapped == 0)); then
    echo -n "no-alignment" 2>&1
    exit 1
fi

error_rate=$(( 100 * mismatches / bases_mapped ))
echo "mismatches=$mismatches, bases_mapped=$bases_mapped, error_rate=$error_rate#"

echo "$mismatches $bases_mapped $error_rate" > {output}

if ((error_rate >= 2 )); then
    touch {wildcards.acc}/error_rate_gt_2pct_NOT_CCS_reads
    echo -n "error-rate-greater-than-2-pct-NOT-CCS-reads" 2>&1
    exit 1
fi
"""

rule call:
    input: bam=rules.bam.output.bam, missmatch=rules.missmatch.output, gap=rules.gaps.output
    output: gvcf="{acc}/{acc}.ref.gatk.vcf"
    log: "LOGS/{acc}.call.log"
    threads: 6
    shell: """

{gatk} HaplotypeCaller -R {ref} -I {input.bam} -O {output.gvcf} \
    --ploidy 1 --minimum-mapping-quality 30 --min-base-quality-score 20
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
{gatk} VariantFiltration \\
    -R {ref} \\
    -V {input} \\
    -O {output} \\
    --filter-name "lowAF" \\
    --filter-expression 'vc.getGenotype("'{wildcards.acc}'").getAD().1.floatValue() / vc.getGenotype("'{wildcards.acc}'").getDP() < 0.15' \\
    --filter-name "lowDP" \\
    --filter-expression 'vc.getGenotype("'{wildcards.acc}'").getDP() < 10' \\
    --filter-name "lowAD" \\
    --filter-expression 'vc.getGenotype("'{wildcards.acc}'").getAD().1 < 5'
"""

rule gatk_pass:
    input: rules.filter_variants.output
    output: "{acc}/{acc}.ref.pass.vcf"
    threads: 1
    log: "LOGS/{acc}.pass.log"
    shell: """
cat {input} | grep ^# > {output}
cat {input} | grep -v ^# | grep PASS | awk '{if ($10 !~ "1/2:") print }' >> {output}
grep -vq "^#" {output} || echo "No snps found"
"""

rule genomecov:
    input: bam=rules.bam.output.bam
    output: "{acc}/{acc}.gatk.bam.avg_cov"
    log: "LOGS/{acc}.genomecov.log"
    shell: """
( bedtools genomecov -d -ibam {input} | awk 'BEGIN {{sum=0}}; {{sum+=$3}}; END{{print sum/NR}}' ) 2>{log} > {output}
"""

rule snpeff:
    input: rules.gatk_pass.output
    output: "{acc}/{acc}.ref.snpeff.vcf"
    log: "LOGS/{acc}.snpeff.log"
    threads: 1
    shell: """
{snpEff} ann \\
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
{snpsift} \\
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
    shell: """
    {vcf_validator} -i {input.snpvcf} 2>{log}
"""
