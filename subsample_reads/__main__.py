#!/usr/bin/env python

from subsample_reads import Plotter, Loader, Mapper, Comparator
from logging import getLogger, basicConfig, INFO
from datetime import datetime as dt
import argparse
import sys
import os

# Ensure log directory exists
os.makedirs("log", exist_ok=True)

# Configure logging
basicConfig(
    filename=dt.now().strftime("log/%H%M%S_%d%m%Y.txt"),
    level=INFO,
    format="{asctime} [{levelname}]: {message}",
    style="{",
    datefmt="%H:%M:%S",
)
logger = getLogger("subsample_reads")

def mapper_mode(args):
    """Chart a distribution of reads from given BAM file."""
    Mapper(
        bam_paths=args.in_bam,
        contig=args.contig,
        start=args.start,
        end=args.end,
        interval_length=args.interval_length,
        interval_count=args.interval_count,
        bed_dir=args.bed_dir,
    )

def sample_mode(args):
    """Sample provided BAM file based on regions in BED file."""
    in_bam = Loader(file=args.in_bam)
    in_bam.sample(
        bed_dir=args.bed_dir,
        bed_file=args.bed,
        main_seed=args.seed,
        out_bam=args.out_bam,
    )
    in_bam.close()

def hlala_mode(args):
    """Sample HLALA outputs based on PRG construction data."""
    in_bam = Loader(
        file=os.path.join(args.hlala_dir, "working", args.sampleID, "remapped_with_a.bam")
    )
    in_bam.hlala(
        hlala_dir=args.hlala_dir,
        bed_dir=args.bed_dir,
        bed_file=args.bed,
        main_seed=args.seed,
        out_bam=args.out_bam,
    )
    in_bam.close()

def compare_mode(args):
    """Compare two BAM files to see how many reads overlap with each other."""
    Comparator(bam1_path=args.bam1, bam2_path=args.bam2, out=args.out)

def plotter_mode(args):
    """Chart a distribution of reads from given BAM file."""
    plotter = Plotter(
        in_bam=args.in_bam,
        map_bam=args.map_bam,
        sub_bam=args.sub_bam,
        bed_dir=args.bed_dir,
        bed=args.bed,
        out_plt=args.out_plt,
    )

    plotter.plot()

def main():
    """Main CLI entry point for subsample-reads toolkit."""
    logger.info("Begin log")

    parser = argparse.ArgumentParser(
        prog="python -m subsample-reads",
        description="Toolkit for mapping, sampling, and plotting BAM files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        required=True,
        title="Functions",
        description="Use -h flag with any subcommand to learn usage.",
    )

    # Mapping
    mapper = subparsers.add_parser(
        "map",
        help="Generate read depth distribution(s) from supplied BAM file(s) and write to BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mapper.add_argument("--in-bam", nargs="+", required=True, help="BAM files to sample read counts from.")
    mapper.add_argument("--contig", required=True, help="A valid contig name present in all provided BAM file(s).")
    mapper.add_argument("--start", required=True, help="Start coordinate of the region to map.")
    mapper.add_argument("--end", required=True, help="End coordinate of the region to map.")
    mapper.add_argument("-d", "--bed-dir", default="bed/", help="Directory for output BED files.")
    intervals = mapper.add_mutually_exclusive_group(required=True)
    intervals.add_argument("--interval-length", default=None, help="Length of intervals to generate.")
    intervals.add_argument("--interval-count", default=None, help="Number of intervals to generate.")
    mapper.set_defaults(func=mapper_mode)

    # Sampling
    sample = subparsers.add_parser(
        "sample",
        help="Apply generated read depth distribution(s) from selected BED file(s) to supplied BAM file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sample.add_argument("--in-bam", required=True, help="Target BAM file to subset.")
    bed_selection = sample.add_mutually_exclusive_group()
    bed_selection.add_argument("--bed-dir", default="bed/", help="Directory to fetch BED files from.")
    bed_selection.add_argument("--bed", default=None, help="Specify one BED file to sample from.")
    sample.add_argument("--seed", default=42, type=int, help="Seed for random subsampling.")
    sample.add_argument("--out-bam", default="out.bam", help="Output BAM file.")
    sample.set_defaults(func=sample_mode)

    # HLA-LA specific sampling
    hlala = subparsers.add_parser(
        "hlala",
        help="Apply generated read depth distribution(s) from selected BED file(s) to HLA-LA output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    hlala.add_argument("--hlala-dir", default="HLA-LA/", help="Path to HLA-LA setup directory.")
    hlala.add_argument("--sampleID", required=True, help="HLA-LA processed sample to subsample from.")
    bed_selection = hlala.add_mutually_exclusive_group()
    bed_selection.add_argument("--bed-dir", default="bed/", help="Directory to fetch BED files from.")
    bed_selection.add_argument("--bed", default=None, help="Specify one BED file to sample from.")
    hlala.add_argument("--seed", default=42, type=int, help="Seed for random subsampling.")
    hlala.add_argument("--out-bam", default="out.bam", help="Output BAM file.")
    hlala.set_defaults(func=hlala_mode)

    # Compare
    compare = subparsers.add_parser(
        "compare",
        help="Compare two BAM files to see how many reads overlap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    compare.add_argument("--bam1", required=True, help="Path to first BAM for the comparison.")
    compare.add_argument("--bam2", required=True, help="Path to second BAM for the comparison.")
    compare.add_argument("--out", required=True, help="Output file for cross-mapping info.")
    compare.set_defaults(func=compare_mode)

    # Plotting
    plotter = subparsers.add_parser(
        "plot",
        help="Plot BAM file(s) read depth together with supplied BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plotter.add_argument("--in-bam", required=True, help="Path to input/original BAM file.")
    plotter.add_argument("--map-bam", required=True, help="Path to mapping BAM file.")
    plotter.add_argument("--sub-bam", required=True, help="Path to subsampled BAM file.")
    bed_selection = plotter.add_mutually_exclusive_group()
    bed_selection.add_argument("--bed-dir", default="bed/", help="Directory to fetch a random BED file from.")
    bed_selection.add_argument("--bed", default=None, help="Specific BED file to plot.")
    plotter.add_argument("--out-plt", default="out.png", help="Path for the output plot.")
    plotter.set_defaults(func=plotter_mode)

    args = parser.parse_args()
    logger.info("Main - Accept arguments")

    try:
        args.func(args)
    except Exception as e:
        logger.error(f"Main - Exception encountered. Details:\n{e}", exc_info=True)
        sys.exit(1)

    logger.info("Main - End log\n")

if __name__ == "__main__":
    main()
