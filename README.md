# subsample-reads
Subsample reads is a Python tool that leverages the `pysam` package to subset fractions of reads across multiple defined chromosomal regions in a single chromosome from a BED file.

Other tools that offer similar functionality such as `samtools view` or `GATK DownsampleSam` only process a single chromosomal region at a time, require multiple command line invocations, as well as the splitting and merging BAM files. `subsample-reads` simplifies this process by allowing users to specify a subsampling pattern across an entire chromosome and outputs a single processed BAM file. This subsampling pattern can also be generated from an existing BAM file.

During subsampling, `subsample-reads` takes into account paired reads that lie across region boundaries and drops paired reads when subsampling in adjacent regions if possible.

The tools contains functions to produce the required BED file that contains sampling regions from a BAM file (`map`), to and plot to quickly check the coverage patterns in BAM files.

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

### Prerequisites:

1. Install Python (tested on Python 3.11.4) dependencies and activate environment before running tool

#### `map`
2. BAM file to produce a read sampling distribution from
3. A valid contig as well as start and end coordinates for the region to map
4. Interval size or count to set up when mapping read coverage across region

#### `sample`
2. BAM to sample according to sampling distribution provided
3. BED file specifying sampling distribution (possibly created by `map`)
4. (Optional) Integer seed to use for sampling (for reproducibility)

#### `plot`
2. One or more BAM files to plot in the provided regions
3. BED file specifying sampling regions (possibly created by `map`)
    - Fraction values are ignored during plotting

---

### Input file specifications:

1. BAM file with the following attributes:
    * Contains a single contig (chromosome or any valid contig name.) Can be generated using `samtools view` from a BAM with wider coverage.
    * Sorted and indexed
    * If multiple BAM files are provided for `plot`, all BAMs should contain reads aligned to the contig provided
2. BED file with the following attributes:
    * No header row
    * Columns representing `contig`, `begin`, `end`, and `fraction` columns in a tab-separated format, not unlike conventional BED files. 
        - `contig` is any contig name featured in the input BAM file header ("chrN" and just "N" contig name formats are both handled internally by `sample`)
        - `begin` is the starting chromosomal coordinate for a sampling interval
        - `end` is the ending coordinate for a sampling interval
        - `fraction` is a fraction in the interval [0.0, 1.0] (inclusive) for the region, all `fraction` values in BED file need to sum to approximately 1.0 (within a +- 0.05 margin).

---

### Execution:

Remember to source the Python virtual environment before running the tool:
```{python}
source venv/bin/activate
```

#### `map`
Executing this command with will produce a BED file with the sampling distribution from the BAM file.
```{python}
python -m subsample_reads map -i <input_BAM_file> -c <contig> -s <region_start> -e <region_end> (-w <window_size>/ -n <window_count>) -r <regions_BED_file>
```

#### `sample`
Executing this command will produce an output file with defined regions being subsetted to the specified fraction of the reads in that region.
```{python}
python -m subsample_reads sample -i <input_BAM_file> -r <regions_BED_file> -s <seed> -o <output_BAM_file>
```

#### `plot`
Executing this command with 1 or more BAM files will produce a coverage plot of the regions in the BED file.
```{python}
python -m subsample_reads plot -i <input_BAM_files> -r <regions_BED_file> -o <output_plot_file>
```
