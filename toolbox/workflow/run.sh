#!/bin/bash -eu

progdir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
progname=$0

platforms=":illumina:ont:pacbio:genbank:"
platform=illumina
instrument=PromethION
accession=
list=
conf=
options=
workdir=workdir
codedir=`dirname $progdir`
extra_params=( '' )

usage() {
cat <<EOF
--DESCRIPTION:
    wrapper script to run SARS2 snakemake

--USAGE:
  $progname --platform --accession|--list [--instrument] [--conf] [--workdir] [--help]
    --platform:   platform, choices = [illumina, ont, pacbio, genbank], default = $platform
    --instrument: applicable to ONT platform only, use with single accession option, default = $instrument
    --accession:  single accession, optional, either --accession or --list must be specified
    --list:       file with accession list, optional, either --accession or --list must be specified
                  NOTE: ONT file list is expected to have two columns <acc> <instrument>
    --conf:       optional custom configfile to override default
    --workdir:    optional working directory, default = $workdir
    --help:       to display this help
EOF
    exit 1
}

while (( $# > 0 ))
do
    case "x$1" in
        x--accession  ) accession="$2";  shift ;;
        x--platform   ) platform="$2";   shift ;;
        x--instrument ) instrument="$2"; shift ;;
        x--list       ) list="$2";       shift ;;
        x--conf       ) conf="$2";       shift ;;
        x--workdir    ) workdir="$2";    shift ;;
        x--help       ) usage                  ;;
        *             ) extra_params+=( $1 )   ;;
    esac
    shift
done

mkdir -p $workdir

if ! [[ "$platforms" =~ ":$platform:" ]]; then
  echo "unknown platform: $platform"
  usage
fi
if [[ -z "$accession" ]] && [[ -z "$list" ]]; then
  echo "please supply either single accession (--accession) or accession list file (--list)"
  usage
fi
if [[ -z "$accession" ]] && [[ -z "$list" ]]; then
  echo "please supply either single accession (--accession) or accession list file (--list), but not both"
  usage
fi
if [[ -z "$accession" ]] && [[ "$platform" = "ont" ]] && ! [[ -z "$instrument" ]]; then
  echo "please supply --instrument for ONT platform"
  usage
fi

## config override options

if ! [[ -z "$conf" ]]; then
  if ! [[ -s "$conf" ]]; then
    echo "configuration file has been specified -${conf}- but does not exists or empty"
    usage
  fi
  options="--configfile $conf"
fi

options="$options --config codedir=\"$codedir\" workdir=$workdir "

if [[ -s "$list" ]]; then
  options="$options accs=$list"
fi
if ! [[ -z "$accession" ]]; then
  echo "$accession $instrument $platform" > $workdir/accs.list
  options="$options accs=accs.list"
fi

snakemake -s $progdir/corona_${platform}.smk -j 6 --printshellcmds $options ${extra_params[*]}
