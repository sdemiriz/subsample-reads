#!/usr/bin/env python

from subsample_reads.Plotter import Plotter
from subsample_reads.Loader import Loader
from subsample_reads.Mapper import Mapper
from subsample_reads.Subsampler import Subsampler
from logging import getLogger, basicConfig, info, DEBUG
import argparse

logger = getLogger("root")

basicConfig(
    filename="log.txt",
    level=DEBUG,
    format="{asctime} [{levelname}]: {message}",
    style="{",
    datefmt="%H:%M:%S",
)


def subsample(args):
    """
    Sample provided BAM file based on regions in BED file
    """
    Subsampler(
        sample_bam_paths=args.in_bams,
        map_bam_paths=args.map_bams,
        bed_dir=args.bed_dir,
        bed_files=args.bed_files,
        bed_count=args.bed_count,
        contig=args.contig,
        start=args.start,
        end=args.end,
        interval_length=args.interval_length,
        interval_count=args.interval_count,
        seed=args.seed,
        out_dir=args.out_dir,
    )


if __name__ == "__main__":

    info("Main - Begin log")

    parser = argparse.ArgumentParser(
        prog="python -m subsample_reads",
        description="Toolkit with functions to map distributions from BAM files, sample BAM files according to supplied distributions and plot resulting read depths across regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--in-bams",
        nargs="+",
        required=True,
        help="One or more BAM files to sample read counts from.",
    )
    parser.add_argument(
        "-m",
        "--map-bams",
        nargs="+",
        required=True,
        help="One or more BAM files to sample read counts from.",
    )
    parser.add_argument(
        "-c",
        "--contig",
        required=True,
        help="A valid contig name present in all provided BAM file(s).",
    )
    parser.add_argument(
        "-s",
        "--start",
        required=True,
        help="The start coordinate of the region to map on the supplied contig.",
    )
    parser.add_argument(
        "-e",
        "--end",
        required=True,
        help="The end coordinate of the region to map on the supplied contig.",
    )

    # Interval length handling group
    handle_intervals = parser.add_mutually_exclusive_group(required=True)
    handle_intervals.add_argument(
        "-l",
        "--interval-length",
        default=None,
        help="Exclusive with --interval-count. Lengths of intervals to generate within supplied start-end region. Final interval may end up shorter due to start-end region length.",
    )
    handle_intervals.add_argument(
        "-n",
        "--interval-count",
        default=None,
        help="Exclusive with --interval-length. Number of intervals to generate within supplied start-end region.",
    )

    # BED file usage handling group
    bed_files = parser.add_mutually_exclusive_group(required=True)
    bed_files.add_argument(
        "-b",
        "--bed-files",
        nargs="+",
        help="Exclusive with --bed-count. Specify the exact BED files to use to generate sampling distribution from.",
    )
    bed_files.add_argument(
        "-n",
        "--bed-count",
        default=1,
        help="Exclusive with --bed-files. Specify how many BED files to reference when generating sampling distribution. Has to be less than or equal to number of BED file(s) present.",
    )

    parser.add_argument(
        "-s",
        "--seed",
        default=42,
        help="Integer seed to direct the random subsampling process.",
    )
    parser.add_argument(
        "-d",
        "--bed-dir",
        default="bed/",
        help="Top level directory name for the output BED files to be placed into, created if does not exist. A subdirectory with the name format contig:start-end will be created based on values supplied to accomodate multiple region selections.",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default="out/",
        help="Directory to write BAM files with subsampled reads and plots.",
    )

    args = parser.parse_args()
    info(f"Main - Accept arguments")

    args.func(args)
    info(f"Main - End log\n")
