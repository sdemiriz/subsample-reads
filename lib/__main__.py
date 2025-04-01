#!/usr/bin/env python

from lib import SubsampleBAM, BAMloader
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Tool to regionally downsample reads in BAM files",
    )

    parser.add_argument("-i", "--in-bam")
    parser.add_argument("-o", "--out-bam")
    parser.add_argument("-r", "--regions")
    parser.add_argument("-S", "--seed")

    args = parser.parse_args()

    print(args)
    bam = BAMloader.BAMloader(in_file=args.in_bam)
    [a for a in bam.get_reads(6, 25_000_000, 35_000_000)]

    bam.plot_pileup(contig=6, start=25_000_000, end=35_000_000)
    # test = SubsampleBAM(
    #     in_file=args.i,
    #     out_file=args.o,
    #     coords=args.c,
    #     seed=args.S,
    #     limits=args.l,
    # )
    # test.run()
