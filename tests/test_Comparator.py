import unittest
import tempfile
import os
import pandas as pd
from subsample_reads.Comparator import Comparator


class TestComparator(unittest.TestCase):
    def setUp(self):
        self.bam1 = "tests/HG00157-HLA-sorted.bam"
        self.bam2 = "tests/HG00157-HLA-sorted.bam"  # self-overlap for simplicity
        self.tempfile = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
        self.out = self.tempfile.name
        self.tempfile.close()

    def tearDown(self):
        if os.path.exists(self.out):
            os.remove(self.out)

    def test_comparator_output(self):
        Comparator(bam1_path=self.bam1, bam2_path=self.bam2, out=self.out)
        self.assertTrue(os.path.isfile(self.out))
        df = pd.read_csv(self.out, sep="\t")
        # Check for expected columns
        expected_cols = [
            "query_name",
            "ref_name_1",
            "ref_start_1",
            "ref_end_1",
            "next_ref_name_1",
            "next_ref_start_1",
            "is_unmapped_1",
            "ref_name_2",
            "ref_start_2",
            "ref_end_2",
            "next_ref_name_2",
            "next_ref_start_2",
            "is_unmapped_2",
            "start_diff",
            "end_diff",
        ]
        for col in expected_cols:
            self.assertIn(col, df.columns)
        self.assertGreater(len(df), 0)

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            Comparator(
                bam1_path="tests/DOES_NOT_EXIST.bam", bam2_path=self.bam2, out=self.out
            )
        with self.assertRaises(FileNotFoundError):
            Comparator(
                bam1_path=self.bam1, bam2_path="tests/DOES_NOT_EXIST.bam", out=self.out
            )


if __name__ == "__main__":
    unittest.main()
