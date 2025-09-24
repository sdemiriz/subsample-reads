#!/bin/bash

#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Script to compare the runtime of different downsampling approach from various tools
# Assumes GATK, samtools, and sambamba are installed and available
# Assumes Python is available and subsample-reads hass all dependencies installed in a virtual environment

# This script will run the downsampling approach for each tool and print the runtime to a file
# The chromosomal coordinates can be adjusted to the region of interest

# Default coordinates are a padded range for the Human Leukocyte Antigen (HLA) region
# By default, the region is divided into intervals of 1000 bases

CHR=chr6
START=25000000
END=35000000
INPUT_BAM=HG002.GRCh38.300x.bam
SEED=42

mkdir -p benchmarks/gatk-inputs

# GATK downsampling cannot limit the operation to a specific interval, operating on the entire BAM file
# The input BAM file is thus pre-processed using samtools first, which is not included in benchmarking
echo "Pre-processing BAM files for GATK..."
INTERVAL_BAM=benchmarks/gatk-inputs/single-interval-subset.bam
samtools view $INPUT_BAM $CHR:$START-$END -b -o $INTERVAL_BAM && \
samtools index $INTERVAL_BAM

echo "Benchmarking GATK HighAccuracy..."
GATK_HIGH_ACCURACY=benchmarks/single-interval-downsamplesam-HighAccuracy.log
{
    time gatk DownsampleSam -P 0.1 -R $SEED -S HighAccuracy -I $INTERVAL_BAM -O /dev/null
} 2> $GATK_HIGH_ACCURACY

echo "Benchmarking GATK ConstantMemory..."
GATK_CONSTANT_MEMORY=benchmarks/single-interval-downsamplesam-ConstantMemory.log
{
    time gatk DownsampleSam -P 0.1 -R $SEED -S ConstantMemory -I $INTERVAL_BAM -O /dev/null
} 2> $GATK_CONSTANT_MEMORY

echo "Benchmarking GATK Chained..."
GATK_CHAINED=benchmarks/single-interval-downsamplesam-Chained.log
{
    time gatk DownsampleSam -P 0.1 -R $SEED -S Chained -I $INTERVAL_BAM -O /dev/null
} 2> $GATK_CHAINED

echo "Benchmarking samtools..."
SAMTOOLS=benchmarks/single-interval-samtools.log
{
    time samtools view -s $SEED.1 -b $INPUT_BAM $CHR:$START-$END -b -o /dev/null
} 2> $SAMTOOLS

echo "Benchmarking sambamba..."
SAMBAMBA=benchmarks/single-interval-sambamba.log
{
    time sambamba view -s $SEED.1 -b $INPUT_BAM $CHR:$START-$END -b -o /dev/null
} 2> $SAMBAMBA

echo "Benchmarking GATK DownsampleByDuplicateSet..."
gatk AddOrReplaceReadGroups -I benchmarks/gatk-inputs/$CHR-$START-$END.bam -O benchmarks/gatk-inputs/$CHR-$START-$END-with-read-groups.bam -LB 1 -PL ILLUMINA -PU unit1 -SM HG002
GATK_DOWN_BY_DUP_SET=benchmarks/single-interval-gatk-DownsampleByDuplicateSet.log
{
    time gatk DownsampleByDuplicateSet --fraction-to-keep 0.1 -I benchmarks/gatk-inputs/$CHR-$START-$END-with-read-groups.bam -O /dev/null
} 2> $GATK_DOWN_BY_DUP_SET

# subsample-reads mapping is also done outside of the benchmarking process
# this assumes the virtual environment with all dependencies is already set up
TEMPLATE=benchmarks/single-interval-benchmark.bed
source venv/bin/activate && \
python -m subsample_reads map --in-bam $INPUT_BAM --contig $CHR --start $START --end $END --interval-count 1 --bed $TEMPLATE

echo "Benchmarking subsample-reads..."
SUBSAMPLE_READS=benchmarks/single-interval-subsample-reads.log
{
    time python -m subsample_reads sample --seed $SEED --in-bam $INPUT_BAM --bed $TEMPLATE --out-bam /dev/null
} 2> $SUBSAMPLE_READS

get_time(){
    grep $1 $2 | cut -d' ' -f2
}

OUTPUT=benchmarks/single-interval-benchmark.txt
make_table() {

    echo "tool(+mode)	real	user	sys" > $OUTPUT
    echo "gatk-HighAccuracy	$(get_time real $GATK_HIGH_ACCURACY)	$(get_time user $GATK_HIGH_ACCURACY)	$(get_time sys $GATK_HIGH_ACCURACY)" >> $OUTPUT
    echo "gatk-ConstantMemory	$(get_time real $GATK_CONSTANT_MEMORY)	$(get_time user $GATK_CONSTANT_MEMORY)	$(get_time sys $GATK_CONSTANT_MEMORY)" >> $OUTPUT
    echo "gatk-Chained	$(get_time real $GATK_CHAINED)	$(get_time user $GATK_CHAINED)	$(get_time sys $GATK_CHAINED)" >> $OUTPUT
    echo "gatk-DownsampleByDuplicateSet	$(get_time real $GATK_DOWN_BY_DUP_SET)	$(get_time user $GATK_DOWN_BY_DUP_SET)	$(get_time sys $GATK_DOWN_BY_DUP_SET)" >> $OUTPUT
    echo "samtools	$(get_time real $SAMTOOLS)	$(get_time user $SAMTOOLS)	$(get_time sys $SAMTOOLS)" >> $OUTPUT
    echo "sambamba	$(get_time real $SAMBAMBA)	$(get_time user $SAMBAMBA)	$(get_time sys $SAMBAMBA)" >> $OUTPUT
    echo "subsample-reads	$(get_time real $SUBSAMPLE_READS)	$(get_time user $SUBSAMPLE_READS)	$(get_time sys $SUBSAMPLE_READS)" >> $OUTPUT
}

make_table