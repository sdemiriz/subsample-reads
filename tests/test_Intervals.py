import unittest
from lib.Intervals import Intervals


class TestIntervals(unittest.TestCase):

    def testNonUniqueContigsInBED(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(file="tests/contigs-not-unique.bed")

        self.assertEqual(
            str(context.exception), f"Not all chr values in BED file are the same"
        )

    def testFractionsOutsideInterval01(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(file="tests/fractions-outside-expected-interval.bed")

        self.assertEqual(
            str(context.exception), f"Fraction values not within [0.0. 1.0] interval"
        )

    def testFractionsNotCloseToOne(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(file="tests/not-sums-to-oneish.bed")

        self.assertEqual(
            str(context.exception), f"Fraction values do not sum close to 1.0"
        )

    def testOverlappingIntervals(self):
        with self.assertRaises(AssertionError) as context:
            Intervals(file="tests/overlapping-intervals.bed")

        self.assertEqual(
            str(context.exception), f"BED file contains overlapping intervals"
        )


if __name__ == "__main__":
    unittest.main()
