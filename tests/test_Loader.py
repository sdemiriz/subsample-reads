import unittest

import pysam

from subsample_reads.Loader import Loader


class TestLoader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_bam = "tests/HG00157-100-mapped-reads.bam"
        cls.bam = Loader(bam_path=cls.test_bam)

    @classmethod
    def tearDownClass(cls):
        cls.bam.close()

    def test_loader_initialization(self):
        self.assertIsInstance(self.bam.bam, pysam.AlignmentFile)
        self.assertEqual(self.bam.bam.filename.decode(), self.test_bam)

    def test_loader_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            Loader(bam_path="tests/NONEXISTENT.bam")

    def test_overlap(self):
        self.assertTrue(self.bam.overlap((0, 10), (0, 10)))  # identity
        self.assertFalse(self.bam.overlap((0, 10), (20, 30)))  # no overlap
        self.assertTrue(self.bam.overlap((0, 5), (3, 7)))  # overlap with overhangs
        self.assertTrue(self.bam.overlap((0, 10), (4, 5)))  # one inside other
        self.assertTrue(self.bam.overlap((0, 10), (9, 20)))  # overlap
        self.assertFalse(self.bam.overlap((0, 10), (10, 20)))  # X end Y start overlap
        self.assertFalse(self.bam.overlap((0, 10), (11, 20)))  # sequential

    def test_normalize_contig(self):
        # Should return the same if already correct
        self.assertEqual(self.bam.normalize_contig("chr6"), "chr6")
        # Should convert '6' to 'chr6' if present
        if "6" in self.bam.bam.references:
            self.assertEqual(self.bam.normalize_contig("6"), "chr6")
        # Should raise if contig is not present
        with self.assertRaises(ValueError):
            self.bam.normalize_contig("chrFAKE")

    def test_get_reference_name(self):
        # Should return 'chr6' for the correct reference id
        ref_id = self.bam.bam.get_tid("chr6")
        self.assertEqual(self.bam.get_reference_name(ref_id), "chr6")

    def test_fetch(self):
        reads = list(self.bam.fetch())
        self.assertTrue(all(r.is_mapped for r in reads))
        self.assertGreater(len(reads), 0)


if __name__ == "__main__":
    unittest.main()
