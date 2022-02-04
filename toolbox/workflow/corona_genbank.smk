

ref = config["ref"]
snpeff_config = config["snpeff_config"]
snpEff2tsv = config["snpEff2tsv"]

products = ["fa", "ref.bam", "ref.bam.bai", "ref.mpileup", "vcf",
            "annotated.vcf", "snpEff_summary.html", "snpEff_summary.genes.txt", "tsv"]

accessions_file_path = config["accs"]
with open(accessions_file_path,'r') as f:
    accessions = [line.strip() for line in f.readlines()]

rule all:
    input: expand("{acc}/{acc}.{product}", acc=accessions, product=products)
    shell: """
ls {input}
"""

rule fetch:
    output: "{acc}/{acc}.fa"
    log: "{acc}/LOGS/{acc}.efetch.log"
    shell: """
        if [[ ! -d {wildcards.acc} ]]; then mkdir {wildcards.acc}; fi
        efetch -db nuccore -id {wildcards.acc} -format fasta > {output} 2> {log}
    """

rule align:
    input: rules.fetch.output
    output: bam = "{acc}/{acc}.ref.bam", index = "{acc}/{acc}.ref.bam.bai"
    threads: 8
    log: log1 = "{acc}/LOGS/{acc}.minimap2.log", log2 = "{acc}/LOGS/{acc}.samtools.view.log",
         log3 = "{acc}/LOGS/{acc}.samtools.sort.log", log4 = "{acc}/LOGS/{acc}.samtools.index.log"
    shell: """
        minimap2 -t {threads} -a {ref} {input} 2> {log.log1} \
            | samtools view -Sb -F4 - 2> {log.log2} \
            | samtools sort - > {wildcards.acc}/{wildcards.acc}.ref.bam 2> {log.log3}
        samtools index {wildcards.acc}/{wildcards.acc}.ref.bam 2> {log.log4}
    """

rule pileup:
    input: rules.align.output.bam
    output: mpileup = "{acc}/{acc}.ref.mpileup"
    log: "{acc}/LOGS/{acc}.mpileup.log"
    shell: """
    bcftools mpileup -a INFO/AD -a FORMAT/AD -a FORMAT/DP --indel-size 125 --min-ireads 1 --max-idepth 1000000 \
        --max-depth 1000000 \
        --fasta-ref {ref} {input} > {output} 2> {log}
    """

rule call:
    input: rules.pileup.output.mpileup
    output: vcf = "{acc}/{acc}.vcf"
    log: "{acc}/LOGS/{acc}.call.log"
    shell: """
    bcftools call -Ov --keep-alts --variants-only --multiallelic-caller -P 0 --ploidy 1 \
        -o {wildcards.acc}/{wildcards.acc}.vcf {input} 2> {log}
    """

rule snpeff:
    input: rules.call.output.vcf
    output: "{acc}/{acc}.snpEff_summary.genes.txt", "{acc}/{acc}.snpEff_summary.html",
            annoVCF = "{acc}/{acc}.annotated.vcf"
    log: "{acc}/LOGS/{acc}.snpeff.log"
    shell: """
    snpeff ann -nodownload -canon -formatEff -classic -no-upstream -no-downstream \
    -c {snpeff_config} sars2_mp -v {input} > {output.annoVCF} \
    -s {wildcards.acc}/{wildcards.acc}.snpEff_summary.html 2> {log}
    """

rule to_tsv:
    input: rules.snpeff.output.annoVCF
    output:  tsv = "{acc}/{acc}.tsv"
    log: "{acc}/LOGS/{acc}.tsv.log"
    shell: """ cat {input} | {snpEff2tsv} > {output.tsv} 2> {log}
    """