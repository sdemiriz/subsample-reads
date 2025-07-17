import unittest
import tempfile
import os
from subsample_reads.Plotter import Plotter


class TestPlotter(unittest.TestCase):
    def setUp(self):
        self.bam = "tests/HG00157-HLA-sorted.bam"
        self.bed = "tests/HG00157-HLA-sorted.bed"
        self.bed_dir = "tests"
        self.tempfile = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
        self.out_plt = self.tempfile.name
        self.tempfile.close()

    def tearDown(self):
        if os.path.exists(self.out_plt):
            os.remove(self.out_plt)

    def test_plotter_output(self):
        plotter = Plotter(
            in_bam=self.bam,
            map_bam=self.bam,
            sub_bam=self.bam,
            bed_dir=self.bed_dir,
            bed=self.bed,
            out_plt=self.out_plt,
        )
        plotter.plot()
        self.assertTrue(os.path.isfile(self.out_plt))
        self.assertGreater(os.path.getsize(self.out_plt), 0)

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            Plotter(
                in_bam="tests/DOES_NOT_EXIST.bam",
                map_bam=self.bam,
                sub_bam=self.bam,
                bed_dir=self.bed_dir,
                bed=self.bed,
                out_plt=self.out_plt,
            )
        with self.assertRaises(FileNotFoundError):
            Plotter(
                in_bam=self.bam,
                map_bam=self.bam,
                sub_bam=self.bam,
                bed_dir=self.bed_dir,
                bed="tests/DOES_NOT_EXIST.bed",
                out_plt=self.out_plt,
            )


if __name__ == "__main__":
    unittest.main()
