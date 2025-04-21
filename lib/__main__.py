#!/usr/bin/env python

from lib.Intervals import Intervals
from lib.BAMloader import BAMloader
from lib.BAMplotter import BAMplotter
from lib.BAMcharter import BAMcharter
import argparse, logging
import numpy as np

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
    logging.info(f"Load input BAM file: {args.in_bam}")

    out_bam = BAMloader(file=args.out_bam, template=in_bam.bam)
    logging.info(f"Open output BAM file: {args.out_bam}")

    bed = Intervals(file=args.regions)
    logging.info(f"Load BED file: {args.regions}")

    # , chr_length=bam.get_length(contig=contig))

    logging.info(f"Start sampling")
    for interval in bed.tree:

        logging.info(
            f"Start subsample interval {bed.contig}:{interval.begin}-{interval.end}"
        )

        in_bam.downsample_reads(
            out_bam=out_bam,
            contig=bed.contig,
            start=interval.begin,
            end=interval.end,
            fraction=interval.data,
            seed=args.seed,
        )

        logging.info(
            f"End subsample interval {bed.contig}:{interval.begin}-{interval.end}"
        )

    in_bam.close()
    out_bam.close()


def plot_mode(args):
    """
    Chart a distribution of reads from given BAM file
    """
    BAMplotter(bam_files=args.in_bam, bed_file=args.regions, out=args.output)


if __name__ == "__main__":

    logging.info("Begin log")

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
    sample.set_defaults(func=plot_mode)

    args = parser.parse_args()
    logging.info(f"Accept arguments: {args}")

    args.func(args)

    logging.info(f"End log\n")

    # BAMplotter([args.in_bam, args.out_bam], 6, 29000000, 29200000)
