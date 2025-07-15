from intervaltree import Interval, IntervalTree
from pathlib import Path
from logging import info
import pandas as pd
import numpy as np


class Intervals:

    def __init__(self, bed_dir: str, bed_file: str) -> None:
        """
        Class constructor: read, validate BED and populate IntervalTree
        """
        info(f"Intervals - Initialize Intervals")

        # Read BED file
        self.handle_bed_files(bed_dir=bed_dir, bed_file=bed_file)
        self.get_bed(path=self.bed_file)

        # Get contig, start, and end from BED
        self.get_contig()
        self.get_limits()

        # Build IntervalTree from BED
        self.populate_tree()
        self.validate()

        info(f"Intervals - Complete initialize Intervals")

    def handle_bed_files(self, bed_dir: str, bed_file: str) -> None:
        """
        Decide which BED files to read in based on provided file list or count
        """
        info("Intervals - Handle BED files")

        if bed_file:
            self.bed_file = bed_file
            info(f"Intervals - Received BED file path {self.bed_file}")
        elif bed_dir:
            self.bed_file = np.random.choice(a=list(Path(bed_dir).glob("*.bed")))
            info(f"Intervals - Selected {self.bed_file} random BED file from {bed_dir}")

    def __len__(self):
        return len(self.tree)

    def get_contig(self):
        """
        Gets contig from first BED file (assumes BED files have been validated)
        """
        self.contig = self.bed["contig"][0]
        info(f"Intervals - Set contig value {self.contig}")

    def get_bed(self, path: str) -> None:
        """
        Read BED file from supplied filename
        """
        info(f"Intervals - Read BED from {path}")

        self.bed = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["contig", "start", "end", "read_count", "fraction"],
            dtype={
                "contig": str,
                "start": int,
                "end": int,
                "read_count": int,
                "fraction": float,
            },
        )

    def populate_tree(self) -> None:
        """
        Populate IntervalTree using rows from BED file
        """
        info("Intervals - Populate interval tree using mean read_counts")

        tree = IntervalTree()
        for row in self.bed.itertuples(index=False):
            tree.add(Interval(begin=row[1], end=row[2], data=row[3]))

        self.tree = IntervalTree(sorted(tree))
        info(f"Intervals - Created IntervalTree from {len(self.tree)} intervals")

    def get_limits(self) -> None:
        """
        Return start and end of region from the first BED file (assumes BED files have been validated)
        """
        self.start, self.end = min(self.bed["start"]), max(self.bed["end"])

    def validate(self) -> None:
        """
        Validate BED files when
        """
        info(f"Intervals - Validate last BED file and tree")

        contig = ""
        start, end = None, None

        # Check for different contigs within BED files
        assert (
            len(pd.unique(self.bed["contig"])) == 1
        ), f"Not all contig values in BED file identical"

        # Check for different contigs between BED files
        if contig:
            assert contig == self.bed["contig"][0], "Contig differs from other contigs"
            contig = self.bed["contig"][0]

        # Check all start and end cooridnates match between BED files
        if start and end:
            assert start == min(
                self.bed["start"]
            ), "Start coordinate differs from other start coordinates"
            assert end == max(
                self.bed["end"]
            ), "End coordinate differs from other end coordinates"
            start, end = min(self.bed["start"]), max(self.bed["end"])

        # Check for interval overlaps
        before = self.tree.copy()
        self.tree.split_overlaps()
        after = self.tree
        assert before == after, "BED file contains overlapping intervals"
