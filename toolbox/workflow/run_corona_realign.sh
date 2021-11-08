#!/bin/bash -eu

ref=${1:-NC_045512.2.fa}
accs=${2:-accs}

[[ -f $ref ]]
[[ -f $acc ]]

echo snakemake -s corona_realign.smk -j 6 --printshellcmds --config ref=$ref accs=$accs
