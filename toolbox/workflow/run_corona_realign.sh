#!/bin/bash -eu

dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

config=${CORONA_REALIGN_CONFIG}

if [ -z $config ]; then
  config=${1:-corona_realign.config.yaml}
fi

[[ -f $config ]]

snakemake -s $dir/corona_realign.smk -j 6 --printshellcmds --configfile $config
# --config ref=$ref accs=$accs
