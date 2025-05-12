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
    for file in args.in_bams:
        Mapper(
            bam_filename=file,
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
        intervals=args.regions,
        initial_seed=args.seed,
        out_bam=args.out_bam,
    )

    in_bam.close()


def plotter_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    Plotter(bam_files=args.in_bam, bed_file=args.regions, out=args.output)


if __name__ == "__main__":

    info("Main - Begin log")

    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Regionally downsample reads in BAM files",
    )

    subparsers = parser.add_subparsers(
        title="subcommands", description="valid subcommands", required=True
    )

    # Mapping
    mapper = subparsers.add_parser("map")
    mapper.add_argument("-i", "--in-bams", nargs="+", required=True)

    intervals = mapper.add_mutually_exclusive_group(required=True)
    intervals.add_argument("-l", "--interval-length", default=None)
    intervals.add_argument("-n", "--interval-count", default=None)

    mapper.add_argument("-c", "--contig", required=True)
    mapper.add_argument("-s", "--start", required=True)
    mapper.add_argument("-e", "--end", required=True)

    mapper.add_argument("-d", "--bed_dir", default="bed/")
    mapper.set_defaults(func=mapper_mode)

    # Sampling
    sample = subparsers.add_parser("sample")
    sample.add_argument("-i", "--in-bam", required=True)
    sample.add_argument("-r", "--regions", required=True)
    sample.add_argument("-o", "--out-bam", default="out.bam")
    sample.add_argument("-s", "--seed", default=42)
    sample.set_defaults(func=sample_mode)

    # Plotting
    plotter = subparsers.add_parser("plot")
    plotter.add_argument("-i", "--in-bam", nargs="+", required=True)
    plotter.add_argument("-r", "--regions", required=True)
    plotter.add_argument("-o", "--output", default="out.png")
    plotter.set_defaults(func=plotter_mode)

    args = parser.parse_args()
    info(f"Main - Accept arguments: {args}")

    args.func(args)

    info(f"End log\n")
