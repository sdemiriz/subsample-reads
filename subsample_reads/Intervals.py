from intervaltree import Interval, IntervalTree
from pathlib import Path
from logging import info
import pandas as pd
import random


class Intervals:

    def __init__(self, files: str) -> None:
        """
        Class constructor: read, validate BED and populate IntervalTree
        """
        info(f"Intervals - Initialize Intervals")
        self.files = files

        # Read and validate BED files
        self.get_beds()
        self.validate()

        # Get contig, start, end and various statistics from BEDs
        self.get_stats()
        self.get_contig()
        self.get_limits()

        # Consolidate BEDs into IntervalTree
        self.populate_tree(query="mean")

        self.total_read_counts = [sum(bed["read_count"]) for bed in self.beds]
        self.fractions = [
            bed["read_count"] / sum(bed["read_count"]) for bed in self.beds
        ]

        info(f"Intervals - Complete initialize Intervals")

    def __len__(self):
        return len(self.beds)

    def get_beds(self):
        """
        Read all input files as DataFrames
        """
        info("Intervals - Read in BED file(s)")
        self.beds = [self.read_bed(path=file) for file in self.files]

    def get_contig(self):
        """
        Gets contig from first BED file (assumes BED files have been validated)
        """
        info("Intervals - Set contig value")
        self.contig = self.beds[0]["contig"][0]

    def get_stats(self):
        """
        Calculate statistics from BED file read counts
        """
        info(f"Intervals - Calculate collective statistics for BED files")

        values = pd.DataFrame()
        for name, bed in zip(self.files, self.beds):
            values[Path(name).stem] = bed["read_count"]

        percentiles = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
        stats = values.transpose().describe(percentiles=percentiles).transpose()

        self.stats = pd.concat(
            [self.beds[0][["contig", "start", "end"]], stats, values], axis=1
        )

    def read_bed(self, path: str) -> list[pd.DataFrame]:
        """
        Read BED file from supplied filename
        """
        info(f"Intervals - Read BED {path}")

        return pd.read_csv(
            path,
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

    def populate_tree(self, query: str = "mean") -> None:
        """
        Populate IntervalTree using rows from BED file
        """
        info("Intervals - Populate interval tree using mean read_counts")

        match query:
            case "mean":
                query_index = list(self.stats.columns).index(query)
                iterable = self.stats[query_index]
            case "random":
                random.seed(42)
                iterable = [
                    random.uniform(i[0], i[1])
                    for i in self.stats[["min", "max"]].itertuples(index=False)
                ]
            case _:
                raise Exception("Invalid query value in Intervals.populate_tree()")

        tree = IntervalTree()
        for row_bed, value in zip(self.beds[0].itertuples(index=False), iterable):
            tree.add(Interval(begin=row_bed[1], end=row_bed[2], data=value))

        self.tree = IntervalTree(sorted(tree))

    def get_limits(self) -> tuple[int, int]:
        """
        Return start and end of region from the first BED file (assumes BED files have been validated)
        """
        self.start, self.end = min(self.beds[0]["start"]), max(self.beds[0]["end"])

    def validate(self) -> None:
        """
        Validate BED files when
        """
        info(f"Intervals - Validate last BED file and tree")

        contig = ""
        start, end = None, None
        for bed in self.beds:
            # Check for different contigs within BED files
            assert (
                len(pd.unique(bed["contig"])) == 1
            ), f"Not all contig values in BED file identical"

            # Check for different contigs between BED files
            if contig:
                assert contig == bed["contig"][0], "Contig differs from other contigs"
                contig = bed["contig"][0]

            # Check all start and end cooridnates match between BED files
            if start and end:
                assert start == min(
                    bed["start"]
                ), "Start coordinate differs from other start coordinates"
                assert end == max(
                    bed["end"]
                ), "End coordinate differs from other end coordinates"
                start, end = min(bed["start"]), max(bed["end"])

        # Check for interval overlaps
        before = self.tree.copy()
        self.tree.split_overlaps()
        after = self.tree
        assert before == after, "BED file contains overlapping intervals"
