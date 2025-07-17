import unittest
import os
import pandas as pd
from subsample_reads.Intervals import Intervals


class TestIntervals(unittest.TestCase):
    def test_init_valid_bed(self):
        i = Intervals(bed_dir="tests", bed_file="tests/test-dataframe-dimensions.bed")
        self.assertEqual(i.contig, "chr6")
        self.assertEqual(i.start, 25000000)
        self.assertEqual(i.end, 35000000)
        self.assertEqual(len(i.tree), 10)

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            Intervals(bed_dir="tests", bed_file="tests/DOES_NOT_EXIST.bed")

    def test_non_unique_contigs(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(bed_dir="tests", bed_file="tests/test-contigs-not-unique.bed")
        self.assertIn(
            "Not all contig values in BED file identical", str(context.exception)
        )

    def test_overlapping_intervals(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(bed_dir="tests", bed_file="tests/test-overlapping-intervals.bed")
        self.assertIn("BED file contains overlapping intervals", str(context.exception))

    def test_get_bed_columns(self):
        i = Intervals(bed_dir="tests", bed_file="tests/test-dataframe-dimensions.bed")
        expected_cols = ["contig", "start", "end", "read_count", "fraction"]
        self.assertListEqual(list(i.bed.columns), expected_cols)


if __name__ == "__main__":
    unittest.main()
