#!/bin/bash -eu

dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

config=${CORONA_GENBANK_CONFIG}

if [ -z $config ]; then
  config=${1:-corona_genbank.config.yaml}
fi

[[ -f $config ]]

snakemake -s $dir/corona_genbank.smk -j 6 --printshellcmds --configfile $config
# --config ref=$ref accs=$accs
