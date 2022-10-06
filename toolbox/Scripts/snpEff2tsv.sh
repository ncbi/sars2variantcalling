#!/bin/bash -eu

# author: zaluninvv

grep -v '^#' | while read chr pos id ref alt qual filter info format a; do
  echo "$info" | tr ';' '\n' | sed -rn 's/^EFF=//p' | tr ',' '\n' | grep -vw 'ORF1ab' | while read eff; do
    #eff=$(echo "$info" | tr ';' '\n' | grep -m1 '^EFF')
    kind=$(echo $eff | sed -r  's/^([A-Z_-]+).*/\1/')
    aa=$(echo "$eff" | cut -f 4 -d'|')
    protein=$(echo "$eff" | cut -f 6 -d'|')
    echo -e "$chr\t$pos\t$ref\t$alt\t$kind\t$aa\t$protein"
  done
done
