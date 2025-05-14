from intervaltree import Interval, IntervalTree
from logging import info
import pandas as pd
import math


class Intervals:

    def __init__(self, files: str) -> None:
        """
        Class constructor: read, validate BED and populate IntervalTree
        """
        info(f"Intervals - Initialize Intervals")
        self.files = files

        self.beds, self.trees = [], []
        for file in self.files:
            bed = self.read_bed(file=file)

            self.beds.append(bed)
            self.tree.append(self.populate(bed=bed))

            self.validate_bed(bed=bed)

        info(f"Intervals - Set contig")
        self.contig = self.beds[0]["contig"][0]

        self.total_read_counts = [sum(bed["read_count"]) for bed in self.beds]
        self.fractions = [
            bed["read_count"] / sum(bed["read_count"]) for bed in self.beds
        ]

        info(f"Intervals - Complete initialize Intervals")

    def read_bed(self, file: str) -> list[pd.DataFrame]:
        """
        Read BED file from supplied filename
        """
        info(f"Intervals - Read BED")

        bed = pd.read_csv(
            file,
            sep="\t",
            header=None,
            names=["contig", "start", "end", "read_count"],
            dtype={
                "contig": str,
                "start": int,
                "end": int,
                "read_count": int,
            },
        )

        info(f"Intervals - Complete read BED")
        return bed

    def populate(self, bed: pd.DataFrame) -> None:
        """
        Populate IntervalTree using rows from BED file
        """
        info(f"Intervals - Populate interval tree")

        tree = IntervalTree()
        for row in bed.itertuples():
            assert (
                len(tree.overlap(begin=row[2], end=row[3])) == 0
            ), "BED file contains overlapping intervals"
            tree.add(Interval(begin=row[2], end=row[3], data=row[4]))

        info(f"Intervals - Complete populate interval tree")
        return sorted(tree)

    def get_limits(self) -> tuple[int, int]:
        """
        Return min and max of the region described in BED file
        """
        return (min(self.beds[0]["start"]), max(self.beds[0]["end"]))

    def validate_bed(self, bed: pd.DataFrame) -> None:
        """
        Checks to validate assumptions when reading intervals from BED file
        """
        info(f"Intervals - Validate BED file")

        assert (
            len(pd.unique(bed["contig"])) == 1
        ), f"Not all contig values in BED file are the same"

        abs_tol = 0.05
        assert math.isclose(
            a=sum(bed["fraction"]), b=1.0, abs_tol=abs_tol
        ), f"Fraction values do not sum close to 1.0"

        info(f"Intervals - Complete validate BED file")
