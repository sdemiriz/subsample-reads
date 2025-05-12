from subsample_reads.Loader import Loader
import unittest


class TestLoader(unittest.TestCase):

    test_bam_filename = "tests/HG00157-HLA-sorted.bam"
    bam = Loader(file=test_bam_filename)

    def testContigNotInReferences(self):

        with self.assertRaises(Exception):
            self.bam.get_length(contig="ImpossibleConfigName")

    def testOverlap(self):
        self.assertTrue(self.bam.overlap((0, 10), (0, 10)))  # identity
        self.assertFalse(self.bam.overlap((0, 10), (20, 30)))  # no overlap

        self.assertTrue(self.bam.overlap((0, 5), (3, 7)))  # overlap with overhangs
        self.assertTrue(self.bam.overlap((0, 10), (4, 5)))  # one inside other

        self.assertTrue(self.bam.overlap((0, 10), (9, 20)))  # overlap
        self.assertFalse(self.bam.overlap((0, 10), (10, 20)))  # X end Y start overlap
        self.assertFalse(self.bam.overlap((0, 10), (11, 20)))  # sequential


if __name__ == "__main__":
    unittest.main()
