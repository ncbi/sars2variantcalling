#!/bin/bash -eu

fq=$1
acc=$2
ref=$3
model=$4
work_dir=$5
corona_ref_acc=$6
nanopore_coverage_threshold=$7
threads=$8
logdir=$9

report_error() {
    sz=$(stat -c %s $1)
    if ((sz<10000)); then
        cat $1 1>&2
    else
        head -c 4550 $1 1>&2
        echo "..."      1>&2
        tail -c 4550 $1 1>&2
    fi
}

log=$logdir/$acc.medaka-consensus.log

sum=$(cat $fq | awk '{if (NR % 4 == 2) sum += length ($0) } END {print sum}')
cov=$(( $sum / 29903 ))

echo -n "sum=$sum cov=$cov"

if (( cov < $nanopore_coverage_threshold )); then
    echo -n "coverage-less-than-threshold-$nanopore_coverage_threshold" 1>&2
    exit 1
fi

if [[ ! -d $acc/$work_dir ]]; then mkdir $acc/$work_dir; fi
if ! medaka_consensus -i $fq -d $ref -m $model -t $threads -o $acc/$work_dir &>$log ; then
    report_error $log
    exit 1
fi

sed "s/$corona_ref_acc/$acc/g" $acc/$work_dir/consensus.fasta > $acc/$acc.$work_dir.consensus.fa
mv $acc/$work_dir/consensus.fasta.gaps_in_draft_coords.bed $acc/$acc.$work_dir.gaps.bed
mv $acc/$work_dir/calls_to_draft.bam $acc/$acc.$work_dir.ref.bam
mv $acc/$work_dir/calls_to_draft.bam.bai $acc/$acc.$work_dir.ref.bam.bai

samtools faidx $acc/$acc.$work_dir.consensus.fa
