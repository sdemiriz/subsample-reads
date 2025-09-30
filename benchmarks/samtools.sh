#!/bin/bash

CHR=$1
START=$2
END=$3
INTERVAL_LENGTH=$4
INPUT_BAM=$5
SEED=$6
INPUTS=$7
OUTPUTS=$8

for i in $(seq $START $INTERVAL_LENGTH $((END-INTERVAL_LENGTH))); do
    samtools view -s $SEED.1 $INPUT_BAM $CHR:$i-$((i+INTERVAL_LENGTH)) -b -o $OUTPUTS/samtools.bam
done