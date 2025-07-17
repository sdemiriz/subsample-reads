# subsample-reads Toolkit

A Python toolkit for mapping, subsampling, comparing, and plotting BAM files using BED-defined intervals. Written for bioinformatics  high-throughput sequencing data and compatible with HLA*LA and its PRG-remapping approach.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/sdemiriz/subsample-reads.git
   cd subsample-reads
   ```
2. Install dependencies (ideally in a virtual environment):
   ```bash
   pip install -r requirements.txt
   ```

## Usage:

## Subcommands

### 1. Mapping

Divide a genomic region into intervals and count the number of reads in each interval.
```bash
python -m subsample_reads map [options]

# Required
--in-bam BAM [BAM ...]      :   One or more BAM files to map.
--contig CONTIG             :   Contig name (e.g., chr6) present in all BAM files.
--start START               :   Start coordinate of the region to map.
--end END                   :   End coordinate of the region to map.

# Mutually exclusive
--interval-length LENGTH    :   Specify interval size.
--interval-count COUNT      :   Specify number of intervals.

# Optional
--bed-dir DIR               :   Output directory for BED files (default: bed/).
```

**Example:**
```bash
python -m subsample_reads map \
  --in-bam sample1.bam sample2.bam \
  --contig chr6 --start 25000000 --end 35000000 \
  --interval-count 10
```

### 2. Sampling

Subsample a BAM file according to a BED-defined read depth distribution.
```bash
python -m subsample_reads sample [options]

# Required 
--in-bam BAM        : BAM file to subsample.

# Optional
--bed-dir DIR       : Directory to fetch BED files from (default: `bed/`).
--bed BED           : Specific BED file to use for sampling.
--seed SEED         : Random seed for reproducibility (default: 42).
--out-bam BAM       : Output BAM file (default: `out.bam`).
```

**Example:**
```bash
python -m subsample_reads sample \
  --in-bam sample1.bam \
  --bed bed/sample1.bed \
  --out-bam subsampled.bam
```

---

### 3. HLA*LA-specific Sampling

Subsample HLA*LA output BAMs (`remapped_with_a.bam`) using PRG-mapped reads and a BED file.
```bash
python -m subsample_reads hlala [options]

# Required
--hlala-dir DIR     : Path to HLA*LA setup directory (default: `HLA-LA/`).
--sampleID ID       : HLA*LA processed sample ID.

# Optional
--bed-dir DIR       : Directory to fetch BED files from (default: `bed/`).
--bed BED           : Specific BED file to use for sampling.
--seed SEED         : Random seed (default: 42).
--out-bam BAM       : Output BAM file (default: `out.bam`).
```

**Example:**
```bash
python -m subsample_reads hlala \
  --hlala-dir HLA-LA/ \
  --sampleID HG00157 \
  --bed bed/sample1.bed \
  --seed 02836 \
  --out-bam hlala-subsampled.bam
```

---

### 4. Comparison

Compare two BAM files to determine how many reads overlap between them. Identify and quantify overlapping reads between two BAM files.
```bash
python -m subsample_reads compare [options]

# Required
--bam1 BAM          : First BAM file.
--bam2 BAM          : Second BAM file.
--out FILE          : Output file for overlap information (tab-delimited).
```

**Example:**
```bash
python -m subsample_reads compare \
  --bam1 sample1.bam \
  --bam2 sample2.bam \
  --out overlap.tsv
```

---

### 5. Plotting
Plot read depth and interval coverage from BAM and BED files. Visualize coverage and read distribution across intervals.

```bash
python -m subsample_reads plot [options]

# Required
--in-bam BAM        : Input/original BAM file.
--map-bam BAM       : Mapping BAM file.
--sub-bam BAM       : Subsampled BAM file.

# Optional
--bed-dir DIR       : Directory to fetch a random BED file from (default: `bed/`).
--bed BED           : Specific BED file to plot.
--out-plt PNG       : Output plot file (default: `out.png`).
```

**Example:**
```bash
python -m subsample_reads plot \
  --in-bam sample1.bam --map-bam sample2.bam --sub-bam subsampled.bam \
  --bed bed/sample1.bed --out-plt coverage.png
```

---

## Logging

Log files are written to the `log/` directory with timestamps for each run.
