from intervaltree import Interval, IntervalTree
from logging import info
import pandas as pd
import math


class Intervals:

    def __init__(self, file: str) -> None:
        """
        Class constructor: read, validate BED and populate IntervalTree
        """
        self.file = file
        info(f"Intervals - Initialize Intervals from file {self.file}")

        self.read_bed()
        self.contig = self.bed["contig"][0]
        info(f"Set contig {self.contig}")

        self.populate()
        self.validate_bed()

    def read_bed(self) -> None:
        """
        Read BED file from supplied filename
        """
        info(f"Intervals - Read BED {self.file}")
        self.bed = pd.read_csv(
            self.file,
            sep="\t",
            header=None,
            names=["contig", "start", "end", "fraction"],
            dtype={"contig": str, "start": int, "end": int, "fraction": float},
        )

    def populate(self) -> None:
        """
        Populate IntervalTree using rows from BED file
        """
        info(f"Intervals - Populate interval tree")
        self.tree = IntervalTree()
        for row in self.bed.itertuples():
            assert (
                len(self.tree.overlap(begin=row[2], end=row[3])) == 0
            ), "BED file contains overlapping intervals"
            self.tree.add(Interval(begin=row[2], end=row[3], data=row[4]))

        self.tree = sorted(self.tree)

    def get_limits(self) -> tuple[int, int]:
        """
        Return min and max of the region described in BED file
        """
        return (min(self.bed["start"]), max(self.bed["end"]))

    def validate_bed(self) -> None:
        """
        Checks to validate assumptions when reading intervals from BED file
        """
        assert (
            len(pd.unique(self.bed["contig"])) == 1
        ), f"Not all contig values in BED file are the same"

        assert (self.bed["fraction"] >= 0.0).all() and (
            self.bed["fraction"] <= 1.0
        ).all(), f"Fraction values not within [0.0. 1.0] interval"

        abs_tol = 0.05
        assert math.isclose(
            a=sum(self.bed["fraction"]), b=1.0, abs_tol=abs_tol
        ), f"Fraction values do not sum close to 1.0"

        info(f"Intervals - Validate BED file")
