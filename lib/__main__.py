#!/usr/bin/env python

from lib.Intervals import Intervals
from lib.BAMloader import BAMloader
from lib.BAMplotter import BAMplotter
import argparse, logging
import numpy as np

if __name__ == "__main__":

    logger = logging.getLogger("root")

    logging.basicConfig(
        filename="log.txt",
        level=logging.DEBUG,
        format="{asctime} - [{levelname}]: {message}",
        style="{",
        datefmt="%H:%M:%S",
    )
    logging.info("Begin log")

    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Tool to regionally downsample reads in BAM files",
    )

    parser.add_argument("-i", "--in-bam", required=True)
    parser.add_argument("-r", "--regions", required=True)
    parser.add_argument("-o", "--out-bam", default="out.bam")
    parser.add_argument("-S", "--seed", required=False, default=42)

    args = parser.parse_args()
    logging.info(f"Accept arguments: {args}")

    # in_bam = BAMloader(file=args.in_bam)
    # logging.info(f"Load input BAM file: {args.in_bam}")

    # out_bam = BAMloader(file=args.out_bam, template=in_bam.bam)
    # logging.info(f"Open output BAM file: {args.out_bam}")

    bed = Intervals(file=args.regions)
    logging.info(f"Load BED file: {args.regions}")

    # , chr_length=bam.get_length(contig=contig))

    # logging.info(f"Start subsampling")
    # for interval in bed.tree:

    #     logging.info(
    #         f"Start subsample interval {bed.contig}:{interval.begin}-{interval.end}"
    #     )

    #     in_bam.downsample_reads(
    #         out_bam=out_bam,
    #         contig=bed.contig,
    #         start=interval.begin,
    #         end=interval.end,
    #         fraction=interval.data,
    #         seed=args.seed,
    #     )

    #     logging.info(
    #         f"End subsample interval {bed.contig}:{interval.begin}-{interval.end}"
    #     )

    # in_bam.close()
    # out_bam.close()

    logging.info(f"End log\n")

    # BAMplotter([args.in_bam, args.out_bam], 6, 29000000, 29200000)
