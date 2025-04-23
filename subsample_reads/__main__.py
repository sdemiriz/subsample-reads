#!/usr/bin/env python

from subsample_reads.Intervals import Intervals
from subsample_reads.BAMloader import BAMloader
from subsample_reads.BAMplotter import BAMplotter
from subsample_reads.BAMcharter import BAMcharter
import argparse, logging
from logging import info

logger = logging.getLogger("root")

logging.basicConfig(
    filename="log.txt",
    level=logging.DEBUG,
    format="{asctime} - [{levelname}]: {message}",
    style="{",
    datefmt="%H:%M:%S",
)


def chart_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    BAMcharter(
        bam_file=args.in_bam,
        bed_file=args.regions,
        contig="chr6",
        window_size=args.window_size,
        window_count=args.window_count,
    )


def sample_mode(args):
    """
    Sample provided BAM file based on regions in BED file
    """
    in_bam = BAMloader(file=args.in_bam)
    out_bam = BAMloader(file=args.out_bam, template=in_bam.bam)

    bed = Intervals(file=args.regions)
    info(f"Load BED file: {args.regions}")

    in_bam.run_subsampling(
        contig=bed.contig, tree=bed.tree, initial_seed=args.seed, out_bam=out_bam
    )

    in_bam.close()
    out_bam.close()


def plot_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    BAMplotter(bam_files=args.in_bam, bed_file=args.regions, out=args.output)


if __name__ == "__main__":

    info("Begin log")

    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Regionally downsample reads in BAM files",
    )

    subparsers = parser.add_subparsers(
        title="subcommands", description="valid subcommands", required=True
    )

    chart = subparsers.add_parser("chart")
    chart.add_argument("-i", "--in-bam", required=True)
    chart.add_argument("-r", "--regions", default="out.bed")
    windows = chart.add_mutually_exclusive_group(required=True)
    windows.add_argument("-w", "--window-size", default=None)
    windows.add_argument("-n", "--window-count", default=None)
    chart.set_defaults(func=chart_mode)

    sample = subparsers.add_parser("sample")
    sample.add_argument("-i", "--in-bam", required=True)
    sample.add_argument("-r", "--regions", required=True)
    sample.add_argument("-o", "--out-bam", default="out.bam")
    sample.add_argument("-s", "--seed", required=False, default=42)
    sample.set_defaults(func=sample_mode)

    plot = subparsers.add_parser("plot")
    plot.add_argument("-i", "--in-bam", nargs="+", required=True)
    plot.add_argument("-r", "--regions", required=True)
    plot.add_argument("-o", "--output", default="out.png")
    plot.set_defaults(func=plot_mode)

    args = parser.parse_args()
    info(f"Accept arguments: {args}")

    args.func(args)

    info(f"End log\n")
