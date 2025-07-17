import unittest
import os
import pandas as pd

from subsample_reads.Mapper import Mapper


class TestMapper(unittest.TestCase):
    def setUp(self):
        self.bam = "tests/HG00157-HLA-sorted.bam"
        self.bed = "tests/HG00157-HLA-sorted-test.bed"
        self.bed_dir = "tests"
        self.bed_out = os.path.join(self.bed_dir, "HG00157-HLA-sorted.bed")
        # Remove output file if it exists
        if os.path.exists(self.bed_out):
            os.remove(self.bed_out)

    def tearDown(self):
        if os.path.exists(self.bed_out):
            os.remove(self.bed_out)

    def test_mapper_count(self):
        m = Mapper(
            bam_paths=[self.bam],
            contig="chr6",
            start="25000000",
            end="35000000",
            interval_count="10",
            interval_length=None,
            bed_dir=self.bed_dir,
        )
        self.assertEqual(m.interval_count, 10)
        self.assertIsNone(m.interval_length)
        self.assertTrue(os.path.isfile(self.bed_out))

    def test_mapper_length(self):
        m = Mapper(
            bam_paths=[self.bam],
            contig="chr6",
            start="25000000",
            end="35000000",
            interval_count=None,
            interval_length="1000000",
            bed_dir=self.bed_dir,
        )
        self.assertEqual(m.interval_length, 1000000)
        self.assertIsNone(m.interval_count)
        self.assertTrue(os.path.isfile(self.bed_out))

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            Mapper(
                bam_paths=["tests/DOES_NOT_EXIST.bam"],
                contig="chr6",
                start="25000000",
                end="35000000",
                interval_count="10",
                interval_length=None,
                bed_dir=self.bed_dir,
            )

    def test_no_interval_args(self):
        with self.assertRaises(ValueError):
            Mapper(
                bam_paths=[self.bam],
                contig="chr6",
                start="25000000",
                end="35000000",
                interval_count=None,
                interval_length=None,
                bed_dir=self.bed_dir,
            )

    def test_bed_output_columns(self):
        Mapper(
            bam_paths=[self.bam],
            contig="chr6",
            start="25000000",
            end="35000000",
            interval_count="10",
            interval_length=None,
            bed_dir=self.bed_dir,
        )
        df = pd.read_csv(self.bed_out, sep="\t", header=None)
        self.assertEqual(df.shape[1], 5)


if __name__ == "__main__":
    unittest.main()
