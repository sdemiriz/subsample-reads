import unittest
from subsample_reads.Intervals import Intervals


class TestIntervals(unittest.TestCase):

    def testCorrectColumns(self):
        i = Intervals(files=["tests/test-dataframe-dimensions.bed"])
        for bed in i.beds:
            self.assertListEqual(
                ["contig", "start", "end", "read_count"], list(bed.columns)
            )

    def testCorrectDimensions(self):
        i = Intervals(files=["tests/test-dataframe-dimensions.bed"])
        for bed in i.beds:
            self.assertTupleEqual((10, 4), bed.shape)

    def testDifferentContigsInBED(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(files=["tests/test-contigs-not-unique.bed"])

        self.assertEqual(
            str(context.exception), f"Not all contig values in BED file identical"
        )

    def testOverlappingIntervals(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(files=["tests/test-overlapping-intervals.bed"])

        self.assertEqual(
            str(context.exception), f"BED file contains overlapping intervals"
        )


if __name__ == "__main__":
    unittest.main()
