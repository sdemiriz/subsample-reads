#!/bin/bash

PUBLICATION_DIR=publication

IN_BAM=HG002.GRCh38.300x.bam
MAP_BAM=HG00157.bam
OUT_BAM=$PUBLICATION_DIR/figure-1-subsample-reads.bam
BED=$PUBLICATION_DIR/figure-1.bed
OUT_PLT=$PUBLICATION_DIR/figure-1.png

CHR=chr6
START=29941260
END=29947000
SEED=42

source venv/bin/activate && \
python -m subsample_reads map --in-bam $MAP_BAM --contig $CHR --start $START  --end $END --interval-count 30 --bed $BED && \
python -m subsample_reads sample --seed $SEED --in-bam $IN_BAM --bed $BED --out-bam $OUT_BAM && \
python -m subsample_reads plot --no-det --in-bam $IN_BAM --map-bam $MAP_BAM --out-bam $OUT_BAM --bed $BED --out-plt $OUT_PLT

SAMTOOLS_OUT_BAM=$PUBLICATION_DIR/figure-1-samtools.bam
SAMTOOLS_OUT_PLT=$PUBLICATION_DIR/figure-1-samtools.png

samtools view -b $IN_BAM $CHR:$START-$END -s $SEED.1 -o $SAMTOOLS_OUT_BAM && \
samtools index $SAMTOOLS_OUT_BAM && \
python -m subsample_reads plot --in-bam $IN_BAM --out-bam $SAMTOOLS_OUT_BAM --bed $BED --out-plt $SAMTOOLS_OUT_PLT