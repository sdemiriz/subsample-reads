#!/usr/bin/env python

from lib.BEDloader import BEDloader
from lib.BAMloader import BAMloader
import argparse
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="subsample-reads",
        description="Tool to regionally downsample reads in BAM files",
    )

    parser.add_argument("-i", "--in-bam", required=True)
    parser.add_argument("-r", "--regions", required=True)
    parser.add_argument("-o", "--out-bam", default="out.bam")
    parser.add_argument("-S", "--seed", required=False, default=42)

    args = parser.parse_args()

    np.random.seed(args.seed)

    in_bam = BAMloader(file=args.in_bam)
    out_bam = BAMloader(file=args.out_bam, template=in_bam.bam)
    bed = BEDloader(in_file=args.regions, chr_length=300_000)

    # , chr_length=bam.get_length(contig=contig))

    for interval in bed.tree:

        in_bam.downsample_reads(
            out_bam=out_bam,
            contig=bed.contig,
            start=interval.begin,
            end=interval.end,
            fraction=interval.data,
        )

    in_bam.bam.close()
    out_bam.bam.close()
