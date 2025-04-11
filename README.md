# subsample-reads

Subsample reads is a Python tool that leverages the `pysam` package to subset fractions of reads across multiple defined chromosomal regions in a single chromosome using a BED file. 

Other tools that offer similar functionality such as `samtools view` or `GATK DownsampleSam` only process a single chromosomal region, require numerous command line invocations, as well as splitting and merging BAM files. `subsample-reads` reduces this burden on the user by allowing them to specify their desired subsampling pattern across an entire chromosome and outputs a single BAM file.

`subsample-reads` also takes into account paired reads that lie across user-defined region boundaries and drops paired reads when subsampling in adjacent regions if possible.

## Installing dependencies:
1. Create a Python virtual environment:
    
    `python -m venv venv`
1. Activate environment:
    
    `source venv/bin/activate`
1. Run `pip`:

    `pip install -r requirements.txt`
1. (Optional) Disable environment when done:

    `deactivate`

## Running:

### Requirements:
1. Install Python (tested on Python 3.11.4) dependencies and activate environment before running tool
1. Single-chromosome BAM file, preferably indexed (tool will index for you but takes longer)
1. BED file with regions defined on the same contig as in BAM and a fourth column with subsetting fractions
1. (Optional) Integer seed, default: 42
1. (Optional) Output filename, default: `out.bam`

---

### Input file specifications:
1. BAM file must feature at most one contig. This file can be produced using `samtools view` from a BAM with wider coverage.
1. BED file that specifies regions and subsampling fractions. This file should have no header, with columns representing `chr`, `start`, `end`, and `fraction` columns in a tab-separated format, not unlike conventional BED files. 
- `chr` is be any contig name featured in the input BAM file header
- `start` is the integer starting coordinate for a region
- `end` is the integer ending coordinate for a region
- `fraction` is the floating point fraction in the interval [0.0, 1.0], inclusive for the region

Note: In cases of overlapping regions, the `fraction` value of the region that is lower in the BED file is applied.

---

### Execution:

Executing tool with all required arguments will produce an output file with defined regions being subsetted to the specified fraction of the reads in that region.

```{python}
source venv/bin/activate # activate virtual environment
python -m lib -i <input_BAM_file> -r <regions_BED_file> -s <seed> -o <output_BAM_file> # run tool
```