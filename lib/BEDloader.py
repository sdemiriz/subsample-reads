import pandas as pd
from intervaltree import Interval, IntervalTree


class BEDloader:

    def __init__(self, in_file: str, chr_length: int):

        self.in_file = in_file
        self.load_bed()

        self.tree = IntervalTree()
        self.populate(chr_length=chr_length)

    def load_bed(self):

        self.bed = pd.read_csv(
            self.in_file,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "fraction"],
            dtype={"chr": int, "start": int, "end": int, "fraction": float},
        )

        assert (
            len(pd.unique(self.bed["chr"])) == 1
        ), f"Not all chr values in BED file are the same"

        self.contig = self.bed["chr"][0]

        assert all(
            self.bed["fraction"] >= 0.0 and self.bed["fraction"] <= 1.0
        ), f"Fraction values not within [0.0. 1.0] interval"

    def populate(self, chr_length: int):

        self.tree.addi(begin=0, end=chr_length, data=1.0)

        for row in self.bed.itertuples():
            self.tree.chop(begin=row[2], end=row[3])
            self.tree.addi(begin=row[2], end=row[3], data=row[4])

        self.tree = sorted(self.tree)
