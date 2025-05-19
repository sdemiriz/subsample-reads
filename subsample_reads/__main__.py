#!/usr/bin/env python

from subsample_reads.Plotter import Plotter
from subsample_reads.Loader import Loader
from subsample_reads.Mapper import Mapper
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


def sample_mode(args):
    """
    Sample provided BAM file based on regions in BED file
    """
    in_bam = Loader(file=args.in_bam)

    in_bam.sample(
        bed_dir=args.bed_dir,
        bed_files=args.bed_files,
        bed_count=args.bed_count,
        initial_seed=args.seed,
        out_bam=args.out_bam,
    )

    in_bam.close()


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
        help="One or more BAM files to sample read counts from. Separate BED files produced for each BAM file, using the same filename with a .bed extenstion.",
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
    parser.add_argument(
        "-d",
        "--bed-dir",
        default="bed/",
        help="Top level directory name for the output BED files to be placed into, created if does not exist. A subdirectory with the name format contig:start-end will be created based on values supplied to accomodate multiple region selections.",
    )
    parser.add_argument(
        "-s",
        "--seed",
        default=42,
        help="Integer seed to direct the random subsampling process.",
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
        "-o",
        "--out-bam",
        default="out.bam",
        help="BAM file to write subsampled reads to.",
    )

    parser.add_argument(
        "-p",
        "--plot-dir",
        default="img/",
        help="Directory path to write resulting plot to.",
    )

    args = parser.parse_args()
    info(f"Main - Accept arguments")

    args.func(args)
    info(f"Main - End log\n")
