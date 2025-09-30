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
    gatk DownsampleSam -P 0.1 -R $SEED -S Chained -I $INPUTS/$CHR-$i-$((i+INTERVAL_LENGTH)).bam -O $OUTPUTS/gatk-chained.bam
done