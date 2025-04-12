import pandas as pd
from intervaltree import Interval, IntervalTree
import logging


class BEDloader:

    def __init__(self, file: str, chr_length: int):
        """
        Constructor from file and optional chr_length
        """
        self.file = file
        logging.info(f"Initialize BEDloader using file {self.file}")
        self.load_bed()

        self.tree = IntervalTree()
        self.populate(chr_length=chr_length)

    def load_bed(self):

        self.bed = pd.read_csv(
            self.file,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "fraction"],
            dtype={"chr": int, "start": int, "end": int, "fraction": float},
        )

        assert (
            len(pd.unique(self.bed["chr"])) == 1
        ), f"Not all chr values in BED file are the same"

        self.contig = self.bed["chr"][0]

        assert (self.bed["fraction"] >= 0.0).all() and (
            self.bed["fraction"] <= 1.0
        ).all(), f"Fraction values not within [0.0. 1.0] interval"

        logging.info(f"Loaded BED {self.file}, contig {self.contig}")

    def populate(self, chr_length: int):

        self.tree.addi(begin=0, end=chr_length, data=1.0)
        logging.info(f"Initialize IntervalTree {self.contig}:0-{chr_length}")

        for row in self.bed.itertuples():
            logging.info(
                f"Add Interval {self.contig}-{row[2]}:{row[3]}, fraction {row[4]}"
            )
            self.tree.chop(begin=row[2], end=row[3])
            self.tree.addi(begin=row[2], end=row[3], data=row[4])

        self.tree = sorted(self.tree)
