#!/usr/bin/env python

from subsample_reads.Intervals import Intervals
from subsample_reads.BAMloader import BAMloader
from subsample_reads.Plotter import Plotter
from subsample_reads.Mapper import Mapper
from logging import getLogger, basicConfig, info, DEBUG
import argparse

logger = getLogger("root")

basicConfig(
    filename="log.txt",
    level=DEBUG,
    format="{asctime} - [{levelname}]: {message}",
    style="{",
    datefmt="%H:%M:%S",
)


def mapper_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    Mapper(
        bam_filename=args.in_bam,
        contig=args.contig,
        start=args.start,
        end=args.end,
        interval_length=args.interval_length,
        interval_count=args.interval_count,
        bed_filename=args.regions,
    )


def sample_mode(args):
    """
    Sample provided BAM file based on regions in BED file
    """
    in_bam = BAMloader(file=args.in_bam)
    out_bam = BAMloader(file=args.out_bam, template=in_bam.bam)

    bed = Intervals(file=args.regions)
    info(f"Load BED file: {args.regions}")

    in_bam.run_sampling(
        contig=bed.contig, tree=bed.tree, initial_seed=args.seed, out_bam=out_bam
    )

    in_bam.close()
    out_bam.close()


def plotter_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    Plotter(bam_files=args.in_bam, bed_file=args.regions, out=args.output)


if __name__ == "__main__":

    info("Begin log")

    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Regionally downsample reads in BAM files",
    )

    subparsers = parser.add_subparsers(
        title="subcommands", description="valid subcommands", required=True
    )

    mapper = subparsers.add_parser("map")
    mapper.add_argument("-i", "--in-bam", required=True)

    intervals = mapper.add_mutually_exclusive_group(required=True)
    intervals.add_argument("-l", "--interval-length", default=None)
    intervals.add_argument("-n", "--interval-count", default=None)

    mapper.add_argument("-c", "--contig", default="chr6")
    mapper.add_argument("-s", "--start", default=25_000_000)
    mapper.add_argument("-e", "--end", default=35_000_000)

    mapper.add_argument("-r", "--regions", default="out.bed")
    mapper.set_defaults(func=mapper_mode)

    sample = subparsers.add_parser("sample")
    sample.add_argument("-i", "--in-bam", required=True)
    sample.add_argument("-r", "--regions", required=True)
    sample.add_argument("-o", "--out-bam", default="out.bam")
    sample.add_argument("-s", "--seed", default=42)
    sample.set_defaults(func=sample_mode)

    plotter = subparsers.add_parser("plot")
    plotter.add_argument("-i", "--in-bam", nargs="+", required=True)
    plotter.add_argument("-r", "--regions", required=True)
    plotter.add_argument("-o", "--output", default="out.png")
    plotter.set_defaults(func=plotter_mode)

    args = parser.parse_args()
    info(f"Accept arguments: {args}")

    args.func(args)

    info(f"End log\n")
