#!/usr/bin/env python

from lib.BEDloader import BEDloader
from lib.BAMloader import BAMloader
import argparse, math
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Tool to regionally downsample reads in BAM files",
    )

    parser.add_argument("-i", "--in-bam", required=True)
    parser.add_argument("-r", "--regions", required=True)
    parser.add_argument("-o", "--out-bam")
    parser.add_argument("-S", "--seed", required=False, default=42)

    args = parser.parse_args()

    contig = 6
    np.random.seed(args.seed)

    bam = BAMloader(in_file=args.in_bam)
    bed = BEDloader(
        in_file=args.regions, chr_length=300_000
    )  # , chr_length=bam.get_length(contig=contig))

    for interval in bed.tree:

        reads1, reads2 = bam.get_read_pairs(
            contig=contig, start=interval.begin, end=interval.end
        )

        for read1, read2 in zip(reads1, reads2):
            pass

        # size = math.ceil(interval.data * sum(1 for r in reads))
        # print(interval, sum(1 for r in reads))

    # bam.plot_pileup(contig=6, start=25_000_000, end=35_000_000)
    # test = SubsampleBAM(
    #     in_file=args.i,
    #     out_file=args.o,
    #     coords=args.c,
    #     seed=args.S,
    #     limits=args.l,
    # )
    # test.run()
