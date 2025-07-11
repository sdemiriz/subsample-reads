#!/usr/bin/env python

from subsample_reads.Plotter import Plotter
from subsample_reads.Loader import Loader
from subsample_reads.Mapper import Mapper
from logging import getLogger, basicConfig, info, DEBUG
from datetime import datetime as dt
import argparse

logger = getLogger("root")

basicConfig(
    filename=dt.now().strftime("log/%H%M%S_%d%m%Y.txt"),
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
        bam_paths=args.in_bam,
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
        bed_file=args.bed,
        main_seed=args.seed,
        out_bam=args.out_bam,
    )

    in_bam.close()


def hlala_mode(args):
    """
    Sample HLALA outputs based on PRG construction data
    """
    in_bam = Loader(
        file=args.hlala_dir + "working/" + args.sampleID + "/remapped_with_a.bam"
    )

    in_bam.hlala(
        hlala_dir=args.hlala_dir,
        bed_dir=args.bed_dir,
        bed_file=args.bed,
        main_seed=args.seed,
        out_bam=args.out_bam,
    )

    in_bam.close()


def plotter_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    Plotter(
        in_bam=args.in_bam,
        map_bam=args.map_bam,
        out_bam=args.out_bam,
        bed_dir=args.bed_dir,
        bed_file=args.bed,
        out_plt=args.out_plt,
    )


def main():

    info("Main - Begin log")

    parser = argparse.ArgumentParser(
        prog="python -m subsample-reads",
        description="Toolkit with functions to map distributions from BAM files, sample BAM files according to supplied distributions and plot resulting read depths across regions.",
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
        help=" Generate read depth distribution(s) from supplied BAM file(s) and write to BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mapper.add_argument(
        "--in-bam",
        nargs="+",
        required=True,
        help="One or more BAM files to sample read counts from. Separate BED files produced for each BAM file, using the same filename with a .bed extenstion.",
    )
    mapper.add_argument(
        "--contig",
        required=True,
        help="A valid contig name present in all provided BAM file(s).",
    )
    mapper.add_argument(
        "--start",
        required=True,
        help="The start coordinate of the region to map on the supplied contig.",
    )
    mapper.add_argument(
        "--end",
        required=True,
        help="The end coordinate of the region to map on the supplied contig.",
    )
    mapper.add_argument(
        "-d",
        "--bed-dir",
        default="bed/",
        help="Top level directory name for the output BED files to be placed into, created if does not exist. A subdirectory with the name format contig:start-end will be created based on values supplied to accomodate multiple region selections.",
    )

    intervals = mapper.add_mutually_exclusive_group(required=True)
    intervals.add_argument(
        "--interval-length",
        default=None,
        help="Exclusive with --interval-count. Lengths of intervals to generate within supplied start-end region. Final interval may end up shorter due to start-end region length.",
    )
    intervals.add_argument(
        "--interval-count",
        default=None,
        help="Exclusive with --interval-length. Number of intervals to generate within supplied start-end region.",
    )
    mapper.set_defaults(func=mapper_mode)

    # Sampling
    sample = subparsers.add_parser(
        "sample",
        help="Apply generated read depth distribution(s) from selected BED file(s) to supplied BAM file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sample.add_argument(
        "--in-bam",
        required=True,
        help="Target BAM file to subset based on the distributions in the selected BED file(s).",
    )

    bed_selection = sample.add_mutually_exclusive_group()
    bed_selection.add_argument(
        "--bed-dir",
        default="bed/",
        help="Top level directory to fetch BED files from. The subdirectory to read from with the name format contig:start-end are determined automatically based on the selected BED file(s).",
    )
    bed_selection.add_argument(
        "--bed",
        default=None,
        help="Specify one BED file to sample from.",
    )

    sample.add_argument(
        "--seed",
        default=42,
        help="Integer seed to direct the random subsampling process.",
    )
    sample.add_argument(
        "--out-bam",
        default="out.bam",
        help="BAM file to write subsampled reads to.",
    )
    sample.set_defaults(func=sample_mode)

    hlala = subparsers.add_parser(
        "hlala",
        help="Apply generated read depth distribution(s) from selected BED file(s) to HLA-LA output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    hlala.add_argument(
        "--hlala-dir",
        default="HLA-LA/",
        help="Path to directory where HLA-LA is has been setup.",
    )

    hlala.add_argument("--sampleID", help="HLA-LA processed sample to subsample from")

    bed_selection = hlala.add_mutually_exclusive_group()
    bed_selection.add_argument(
        "--bed-dir",
        default="bed/",
        help="Top level directory to fetch BED files from. The subdirectory to read from with the name format contig:start-end are determined automatically based on the selected BED file(s).",
    )
    bed_selection.add_argument(
        "--bed",
        default=None,
        help="Specify one BED file to sample from.",
    )

    hlala.add_argument(
        "--seed",
        default=42,
        help="Integer seed to direct the random subsampling process.",
    )
    hlala.add_argument(
        "--out-bam",
        default="out.bam",
        help="BAM file to write subsampled reads to.",
    )
    hlala.set_defaults(func=hlala_mode)

    # Plotting
    plotter = subparsers.add_parser(
        "plot",
        help="Plot BAM file(s) read depth together with supplied BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plotter.add_argument(
        "--in-bam",
        required=True,
        help="One or more BAM files to plot depth for based on provided BED file(s).",
    )
    plotter.add_argument(
        "--map-bam",
        required=True,
        help="One or more BAM files to plot depth for based on provided BED file(s).",
    )
    plotter.add_argument(
        "--out-bam",
        required=True,
        help="One or more BAM files to plot depth for based on provided BED file(s).",
    )

    bed_selection = plotter.add_mutually_exclusive_group()
    bed_selection.add_argument(
        "--bed-dir",
        default="bed/",
        help="Top level directory to fetch a BED file from. The subdirectory to read from with the name format contig:start-end are determined automatically based on the selected BED file(s).",
    )
    bed_selection.add_argument(
        "--bed",
        default=None,
        help="TODO",
    )

    plotter.add_argument(
        "--out-plt", default="out.png", help="Path to write resulting plot to."
    )
    plotter.set_defaults(func=plotter_mode)

    args = parser.parse_args()
    info(f"Main - Accept arguments")

    args.func(args)
    info(f"Main - End log\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        info(f"Main - Exception encountered. Details:\n{e}")
        raise e
