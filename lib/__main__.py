#!/usr/bin/env python

from lib import SubsampleBAM
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
    # test = SubsampleBAM(
    #     in_file=args.i,
    #     out_file=args.o,
    #     coords=args.c,
    #     seed=args.S,
    #     limits=args.l,
    # )
    # test.run()
