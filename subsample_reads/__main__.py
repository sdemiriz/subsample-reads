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


def mapper_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    Mapper(
        bam_paths=args.in_bams,
        contig=args.contig,
        start=args.start,
        end=args.end,
        interval_length=args.interval_length,
        interval_count=args.interval_count,
        bed_dir=args.bed_dir,
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


def plotter_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    Plotter(
        bam_files=args.in_bam,
        bed_dir=args.bed_dir,
        bed_files=args.bed_files,
        bed_count=args.bed_count,
        out=args.output,
    )


if __name__ == "__main__":

    info("Main - Begin log")

    parser = argparse.ArgumentParser(
        prog="python -m subsample-reads",
        description="Toolkit with functions to map distributions from BAM files, sample BAM files according to supplied distributions and plot resulting read depths across regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        required=True,
        title="Functions",
        description="Use -h flag with any subcommand to learn more about usage.",
    )

    # Mapping
    mapper = subparsers.add_parser(
        "map",
        help=" Generate read depth distribution(s) from supplied BAM file(s) and write to BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mapper.add_argument(
        "-i",
        "--in-bams",
        nargs="+",
        required=True,
        help="One or more BAM files to sample read counts from. Separate BED files produced for each BAM file, using the same filename with a .bed extenstion.",
    )
    mapper.add_argument(
        "-c",
        "--contig",
        required=True,
        help="A valid contig name present in all provided BAM file(s).",
    )
    mapper.add_argument(
        "-s",
        "--start",
        required=True,
        help="The start coordinate of the region to map on the supplied contig.",
    )
    mapper.add_argument(
        "-e",
        "--end",
        required=True,
        help="The end coordinate of the region to map on the supplied contig.",
    )

    intervals = mapper.add_mutually_exclusive_group(required=True)
    intervals.add_argument(
        "-l",
        "--interval-length",
        default=None,
        help="Exclusive with --interval-count. Lengths of intervals to generate within supplied start-end region. Final interval may end up shorter due to start-end region length.",
    )
    intervals.add_argument(
        "-n",
        "--interval-count",
        default=None,
        help="Exclusive with --interval-length. Number of intervals to generate within supplied start-end region.",
    )
    mapper.add_argument(
        "-d",
        "--bed-dir",
        default="bed/",
        help="Top level directory name for the output BED files to be placed into, created if does not exist. A subdirectory with the name format contig:start-end will be created based on values supplied to accomodate multiple region selections.",
    )
    mapper.set_defaults(func=mapper_mode)

    # Sampling
    sample = subparsers.add_parser(
        "sample",
        help="Apply generated read depth distribution(s) from selected BED file(s) to supplied BAM file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sample.add_argument(
        "-i",
        "--in-bam",
        required=True,
        help="Target BAM file to subset based on the distributions in the selected BED file(s).",
    )
    sample.add_argument(
        "-d",
        "--bed-dir",
        default="bed/",
        help="Top level directory to fetch BED files from. The subdirectory to read from with the name format contig:start-end are determined automatically based on the selected BED file(s).",
    )
    sample.add_argument(
        "-s",
        "--seed",
        default=42,
        help="Integer seed to direct the random subsampling process.",
    )

    bed_files = sample.add_mutually_exclusive_group(required=True)
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

    sample.add_argument(
        "-o",
        "--out-bam",
        default="out.bam",
        help="BAM file to write subsampled reads to.",
    )
    sample.set_defaults(func=sample_mode)

    # Plotting
    plotter = subparsers.add_parser(
        "plot",
        help="Plot BAM file(s) read depth together with supplied BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plotter.add_argument(
        "-i",
        "--in-bam",
        nargs="+",
        required=True,
        help="One or more BAM files to plot depth for based on provided BED file(s).",
    )
    plotter.add_argument(
        "-d",
        "--bed-dir",
        default="bed/",
        help="Top level directory to fetch BED files from. The subdirectory to read from with the name format contig:start-end are determined automatically based on the selected BED file(s).",
    )

    bed_files = plotter.add_mutually_exclusive_group(required=True)
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
        help="Exclusive with --bed-files. Specify how many BED files to reference when generating sampling distribution.",
    )

    plotter.add_argument(
        "-o", "--output", default="out.png", help="Path to write resulting plot to."
    )
    plotter.set_defaults(func=plotter_mode)

    args = parser.parse_args()
    info(f"Main - Accept arguments")

    args.func(args)
    info(f"Main - End log\n")
