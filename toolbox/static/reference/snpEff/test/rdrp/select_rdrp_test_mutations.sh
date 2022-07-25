#!/bin/bash -eu
# extract RdRp mutations around the slippage position from BigQuery as a VCF.

echo '''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=NC_045512.2,length=29903>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##bcftools_normVersion=1.9+htslib-1.9
##bcftools_normCommand=norm --multiallelics -both --fasta-ref NC_045512.2.fa test.vcf; Date=Wed Mar 10 20:06:03 2021
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''

sql="""
select 'NC_045512.2' CHROM, POS, '.' ID, REF, ALT, 255 QUAL, '.' FILTER,  
concat('EFFECT=', coalesce(any_value(EFFECT), ''), ';variation=', coalesce(any_value(variation), ''), ';ref_aa=', coalesce(any_value(ref_aa), ''), ';alt_aa=', coalesce(any_value(alt_aa), ''), ';protein_position=', coalesce(any_value(protein_position), 0)) INFO 
from \`ncbi-sra-virus-variation-prod.SARS_COV_2.annotated_variations\`
where protein_name='RNA-dependent_RNA_polymerase-gene' and POS between 13465 and 13471
group by POS, REF, ALT
order by POS
"""

echo "$sql" | bq --project_id ncbi-sra-virus-variation-prod --quiet query -n=0 --use_legacy_sql=false --format csv --max_rows 100000
