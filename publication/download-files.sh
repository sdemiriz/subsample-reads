#!/bin/bash

echo "Downloading files for publication..."
echo "samtools is required for postprocessing of downloaded files"
echo "Downloading large files, this may take a while"

# Genome In A Bottle WGS 300x sample HG002, download and index
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam && \
samtools index HG002.GRCh38.300x.bam && \

# Thousand Genomes Project sample HG00157, download, convert to BAM and index
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240157/HG00157.final.cram && \
samtools view HG00157.final.cram -b -o HG00157.bam
samtools index HG00157.bam