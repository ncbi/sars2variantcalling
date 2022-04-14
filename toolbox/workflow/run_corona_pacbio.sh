#!/bin/bash -eu

dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

acc=
cpu=6
cfg=corona_pacbio.config.yaml

while (( $# > 0 )); do
  case "x$1" in
    x--cfg ) cfg=$2; shift;;
    x--acc ) acc=$2;    shift;;
    x--cpu ) cpu=$2;    shift;;
  esac
  shift
done

if ! [[ -f $cfg ]] && [[ -f "$dir/$cfg" ]]; then
  cfg=$dir/$cfg
else
  echo "could not find config neither in . nor in $dir"
	usage
fi
if [[ -z "$acc" ]]; then
	echo "ERROR: accession was not supplied, quitting"
	usage
fi

echo $acc>accs
snakemake -s $dir/corona_pacbio.smk --cores $cpu --printshellcmds --configfile $cfg $*
