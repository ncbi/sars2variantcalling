#!/bin/bash -eu
# given a test VCF with snpeff annotations from BigQuery, re-run annotation and print out the mistmacthes.

cat test_rdrp.vcf | snpeff ann -nodownload -formatEff -classic -noStats -noLog -quiet -no-upstream -no-downstream sars2 \
| grep -v '^#' | while read CHROM POS ID REF ALT QUAL FILTER INFO; do 

  variation=$(echo $INFO | sed -r 's/.*variation=([^;]+).*/\1/')
  v2=$(echo $INFO | sed -rn 's/.*EFF=(.*)$/\1/p'  | cut -f 4 -d '|')
  
  if [[ $variation != $v2 ]]; then
    echo -e "Mismatch at $POS:$REF->$ALT: $variation\t$v2"
  fi
done

# example line:
# NC_045512.2	13465	.	A	G	255.0	.	EFFECT=SYNONYMOUS_CODING;variation=L8;ref_aa=L;alt_aa=;protein_position=8;EFF=SYNONYMOUS_CODING(LOW|SILENT|ttA/ttG|L8|930|RNA-dependent_RNA_polymerase-gene|protein_coding|CODING|TRANSCRIPT_RNA-dependent_RNA_polymerase-gene|1|G|WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS)
