# subsample-reads
Subsample reads is a Python tool that leverages the `pysam` package to subset fractions of reads across multiple defined chromosomal intervals in a single contig..

Other tools that offer similar functionality such as `samtools view` or `GATK DownsampleSam` only process a single chromosomal region at a time, require multiple command line invocations and perform downsampling based on only a single fraction. What sets `subsample-reads` apart from these tools is the ability to specify a subsampling distribution across an entire contig based on existing BAM files and apply the generated distribution to other BAM files by sampling a number of reads from each corresponding region to mimic the distribution. The results are visualized separately for each BAM file sampled this way.

During subsampling, `subsample-reads` randomly attributes reads that lie on the boundary of two intervals to either interval. No other checks are currently applied.

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

1. Install Python (tested on 3.11.4) dependencies and activate environment
1. One or more BAM files (sorted, indexed) to map a sampling distribution from
1. A valid contig, start and end coordinates for the region to map
1. Interval lengths or count to define interval sizes across the mapping region
1. One or more BAM files (sorted, indexed) to sample reads based on sampling distribution
1. (Optional) Directories to place mapping data and results into
1. (Optional) Integer seed to use for sampling (for reproducibility)

---

### Execution:

Remember to source the Python virtual environment before running the tool:
```{python}
source venv/bin/activate
python -m subsample_reads \
    -i <BAMs_to_sample_from> \
    -m <BAMs_used_for_mapping> \
    -d <BED_directory> \
    -b <BED_files_to_use_for_sampling> \
    -z <BED_file_count_to_randomly_sample_from> \
    -c <contig_coordinate> \
    -s <start_coordinate> \
    -e <end_coordinate> \
    -i <interval_length> \
    -i <interval_count> \
    -s <seed> \
    -o <output_dir_for_BAMs_and_plots>
```

