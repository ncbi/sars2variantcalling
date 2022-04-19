#!/bin/bash -eu

dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
progname=`basename $0`

acc=
cpu=6
cfg=corona_pacbio.config.yaml

function usage {
cat <<EOF
--USAGE:
  $dir/$progname --acc <accession> [--cfg <config yaml>] [--cpu <number of cores>]
  where:
    --acc <acession>        : required
    --cpu <number of cores> : optional, default = $cpu
    --cfg <path to config>  : optional, default = $cfg
EOF
    exit 1
}

extra_params=( '' )

while (( $# > 0 )); do
  case "x$1" in
    x--cfg ) cfg=$2; shift;;
    x--acc ) acc=$2; shift;;
    x--cpu ) cpu=$2; shift;;
    *      ) extra_params+=( $1 ) ;;
  esac
  shift;
done

#echo ${extra_params[*]}

if ! [[ -f $cfg ]] && [[ -f "$dir/$cfg" ]]; then
  cfg=$dir/$cfg
else
  echo "could not find config $cfg neither in . nor in $dir"
	usage
fi
if [[ -z "$acc" ]]; then
	echo "ERROR: accession was not supplied, quitting"
	usage
fi

echo $acc>accs
snakemake -s $dir/corona_pacbio.smk --cores $cpu --printshellcmds --configfile $cfg ${extra_params[*]}
