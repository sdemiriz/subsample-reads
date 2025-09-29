#!/bin/bash

#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

# Compare the runtime of different downsampling approach from various tools
# Assumes GATK, samtools, and sambamba are installed and available
# Assumes Python is available and subsample-reads hass all dependencies installed in a virtual environment

# This script will run the downsampling approach for each tool and print the runtime to a file
# The chromosomal coordinates can be adjusted to the region of interest

# Default coordinates are a padded range for the Human Leukocyte Antigen (HLA) region
# By default, the region is divided into intervals of 1000 bases

CHR=chr6
START=25000000
END=35000000
INTERVAL_LENGTH=1000
INPUT_BAM=HG002.GRCh38.300x.bam
SEED=42

INPUTS=benchmarks/gatk-inputs
OUTPUTS=benchmarks/outputs
mkdir -p $INPUTS
mkdir -p $OUTPUTS

# GATK downsampling cannot limit the operation to a specific interval, operating on the entire BAM file
# The input BAM file is thus pre-processed using samtools first, which is not included in benchmarking
echo "Pre-processing BAM files for GATK..."
for i in $(seq ${START} ${INTERVAL_LENGTH} $((END-INTERVAL_LENGTH)));
do
    INTERVAL_BAM=$INPUTS/${CHR}-${i}-$((i+INTERVAL_LENGTH)).bam
    samtools view $INPUT_BAM $CHR:$i-$((i+INTERVAL_LENGTH)) -b -o $INTERVAL_BAM && \
    samtools index $INTERVAL_BAM
done

echo "Benchmarking GATK HighAccuracy..."
GATK_HIGH_ACCURACY=$OUTPUTS/multi-interval-gatk-high-accuracy.log
env time -o $GATK_HIGH_ACCURACY -v bash -c "
    for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
        gatk DownsampleSam -P 0.1 -R \$SEED -S HighAccuracy -I benchmarks/gatk-inputs/\${CHR}-\${i}-\$((i+INTERVAL_LENGTH)).bam -O /benchmarks/outputs/gatk-high-accuracy.bam
    done
"

echo "Benchmarking GATK ConstantMemory..."
GATK_CONSTANT_MEMORY=$OUTPUTS/multi-interval-gatk-constant-memory.log
env time -o $GATK_CONSTANT_MEMORY -v bash -c "
    for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
        gatk DownsampleSam -P 0.1 -R \$SEED -S ConstantMemory -I benchmarks/gatk-inputs/\${CHR}-\${i}-\$((i+INTERVAL_LENGTH)).bam -O /benchmarks/outputs/gatk-constant-memory.bam
    done
"

echo "Benchmarking GATK Chained..."
GATK_CHAINED=$OUTPUTS/multi-interval-gatk-chained.log
env time -o $GATK_CHAINED -v bash -c "
        for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
            gatk DownsampleSam -P 0.1 -R \$SEED -S Chained -I benchmarks/gatk-inputs/\${CHR}-\${i}-\$((i+INTERVAL_LENGTH)).bam -O /benchmarks/outputs/gatk-chained.bam
        done
    "

for i in $(seq ${START} ${INTERVAL_LENGTH} $((END-INTERVAL_LENGTH)));
do
    gatk AddOrReplaceReadGroups -I $INPUTS/$CHR-$i-$((i+INTERVAL_LENGTH)).bam -O $INPUTS/$CHR-$i-$((i+INTERVAL_LENGTH))-with-read-groups.bam -LB 1 -PL ILLUMINA -PU unit1 -SM HG002
done

# echo "Benchmarking GATK DownsampleByDuplicateSet..."
# GATK_DOWN_BY_DUP_SET=$OUTPUTS/multi-interval-gatk-downsample-by-dup-set.log
# env time -o $GATK_DOWN_BY_DUP_SET -v bash -c "
#         for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
#             gatk DownsampleByDuplicateSet --fraction-to-keep 0.1 -I $INPUTS/\$CHR-\$i-\$((i+INTERVAL_LENGTH))-with-read-groups.bam -O $OUTPUTS/gatk-downsample-by-dup-set.bam
#         done
#     "

echo "Benchmarking samtools..."
SAMTOOLS=$OUTPUTS/multi-interval-samtools.log
env time -o $SAMTOOLS -v bash -c "
        for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
            samtools view -s \$SEED.1 -b \$INPUT_BAM \$CHR:\$i-\$((i+INTERVAL_LENGTH)) -b -o /benchmarks/outputs/samtools.bam
        done
    "

SAMBAMBA=$OUTPUTS/multi-interval-sambamba.log
echo "Benchmarking sambamba..."
env time -o $SAMBAMBA -v bash -c "
        for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
            sambamba view -s 42.1 \$INPUT_BAM \$CHR:\$i-\$((i+INTERVAL_LENGTH)) -o /benchmarks/outputs/sambamba.bam
        done
    "

# subsample-reads mapping is also done outside of the benchmarking process
# this assumes the virtual environment with all dependencies is already set up
source venv/bin/activate && \
python -m subsample_reads map --in-bam $INPUT_BAM --contig $CHR --start $START --end $END --interval-length $INTERVAL_LENGTH --bed $OUTPUTS/multi-interval-benchmark.bed && \

echo "Benchmarking subsample-reads..."
SUBSAMPLE_READS=$OUTPUTS/multi-interval-subsample-reads.log
env time -o $SUBSAMPLE_READS -v bash -c "
    python -m subsample_reads sample --seed $SEED --in-bam $INPUT_BAM--bed $OUTPUTS/multi-interval-benchmark.bed --out-bam $OUTPUTS/subsample-reads.bam
"

get_time(){
    grep $1 $2 | rev | cut -d' ' -f1 | rev
}

OUTPUT=benchmarks/multi-interval-benchmark.txt
make_table() {

    echo "Tool/Mode	User time	System time	Wall clock  Memory used" > $OUTPUT
    echo "gatk-HighAccuracy	$(get_stat User $GATK_HIGH_ACCURACY)	$(get_stat System $GATK_HIGH_ACCURACY)	$(get_stat wall $GATK_HIGH_ACCURACY)    $(get_stat Maximum $GATK_HIGH_ACCURACY)KB" >> $OUTPUT
    echo "gatk-ConstantMemory	$(get_stat User $GATK_CONSTANT_MEMORY)	$(get_stat System $GATK_CONSTANT_MEMORY)	$(get_stat wall $GATK_CONSTANT_MEMORY)  $(get_stat Maximum $GATK_CONSTANT_MEMORY)KB" >> $OUTPUT
    echo "gatk-Chained	$(get_stat User $GATK_CHAINED)	$(get_stat System $GATK_CHAINED)	$(get_stat wall $GATK_CHAINED)  $(get_stat Maximum $GATK_CHAINED)KB" >> $OUTPUT
    # echo "gatk-DownsampleByDuplicateSet	$(get_stat User $GATK_DOWN_BY_DUP_SET)	$(get_stat System $GATK_DOWN_BY_DUP_SET)	$(get_stat wall $GATK_DOWN_BY_DUP_SET)  $(get_stat Maximum $GATK_DOWN_BY_DUP_SET)KB" >> $OUTPUT
    echo "samtools	$(get_stat User $SAMTOOLS)	$(get_stat System $SAMTOOLS)	$(get_stat wall $SAMTOOLS   $(get_stat Maximum $SAMTOOLS))KB" >> $OUTPUT
    echo "sambamba	$(get_stat User $SAMBAMBA)	$(get_stat System $SAMBAMBA)	$(get_stat wall $SAMBAMBA)  $(get_stat Maximum $SAMBAMBA)KB" >> $OUTPUT
    echo "subsample-reads	$(get_stat User $SUBSAMPLE_READS)	$(get_stat System $SUBSAMPLE_READS)	$(get_stat wall $SUBSAMPLE_READS) $(get_stat Maximum $SUBSAMPLE_READS)KB" >> $OUTPUT
}

make_table