#!/bin/bash -eux

dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
snpeff_home=${dir}/..
data=${snpeff_home}/data
db=sars2
if [[ ! -d $data/$db ]]; then mkdir -p ${data}/${db} ; fi

cp GCF_009858895.2_ASM985889v3_genomic_mature_peptides.gff ${data}/${db}/genes.gff
gzip -f ${data}/${db}/genes.gff
snpEff build -config ${snpeff_home}/snpEff.config -gff3 -v ${db}
cat test.norm.vcf | snpEff ann -nodownload -config ${snpeff_home}/snpEff.config -canon -formatEff -no-upstream -classic -no-downstream sars2 > test.norm.vcf.snpEff
