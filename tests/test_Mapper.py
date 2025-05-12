from subsample_reads.Mapper import Mapper
import unittest, os
import pandas as pd


class TestMapper(unittest.TestCase):

    filename_root = "tests/HG00157-HLA-sorted"
    contig = "chr6"
    start = "25000000"
    end = "35000000"
    count = "10"
    length = "1000000"

    count_mapper = Mapper(
        bam_filenames=[filename_root + ".bam"],
        contig=contig,
        start=start,
        end=end,
        interval_count=count,
        interval_length=None,
        bed_filename=filename_root + ".bed",
    )

    def testCorrectArgumentTypes_Count(self):

        self.assertIsInstance(self.count_mapper.bam_filenames, list)
        self.assertIsInstance(self.count_mapper.contig, str)
        self.assertIsInstance(self.count_mapper.start, int)
        self.assertIsInstance(self.count_mapper.end, int)
        self.assertIsInstance(self.count_mapper.interval_count, int)
        self.assertIsInstance(self.count_mapper.bed_filename, str)

        self.assertIsNone(self.count_mapper.interval_length)

    def testCorrectArgumentValues_Count(self):

        self.assertEqual(self.count_mapper.bam_filenames, [self.filename_root + ".bam"])
        self.assertEqual(self.count_mapper.contig, "chr6")
        self.assertEqual(self.count_mapper.start, int(self.start))
        self.assertEqual(self.count_mapper.end, int(self.end))
        self.assertEqual(self.count_mapper.interval_count, int(self.count))
        self.assertEqual(self.count_mapper.bed_filename, self.filename_root + ".bed")

    length_mapper = Mapper(
        bam_filenames=[filename_root + ".bam"],
        contig=contig,
        start=start,
        end=end,
        interval_count=None,
        interval_length=length,
        bed_filename=filename_root + ".bed",
    )

    def testCorrectArgumentTypes_Length(self):

        self.assertIsInstance(self.length_mapper.bam_filenames, list)
        self.assertIsInstance(self.length_mapper.contig, str)
        self.assertIsInstance(self.length_mapper.start, int)
        self.assertIsInstance(self.length_mapper.end, int)
        self.assertIsInstance(self.length_mapper.interval_length, int)
        self.assertIsInstance(self.length_mapper.bed_filename, str)

        self.assertIsNone(self.length_mapper.interval_count)

    def testCorrectArgumentValues_Length(self):

        self.assertEqual(
            self.length_mapper.bam_filenames, [self.filename_root + ".bam"]
        )
        self.assertEqual(self.length_mapper.contig, "chr6")
        self.assertEqual(self.length_mapper.start, int(self.start))
        self.assertEqual(self.length_mapper.end, int(self.end))
        self.assertEqual(self.length_mapper.interval_length, int(self.length))
        self.assertEqual(self.length_mapper.bed_filename, self.filename_root + ".bed")

    def testBEDcreated(self):

        os.remove(self.filename_root + ".bed")

        Mapper(
            bam_filenames=[self.filename_root + ".bam"],
            contig=self.contig,
            start=self.start,
            end=self.end,
            interval_count=None,
            interval_length=self.length,
            bed_filename=self.filename_root + ".bed",
        )

        self.assertTrue(os.path.isfile(self.filename_root + ".bed"))

    def testBEDcontentsIdentical(self):

        Mapper(
            bam_filenames=[self.filename_root + ".bam"],
            contig=self.contig,
            start=self.start,
            end=self.end,
            interval_count=None,
            interval_length=self.length,
            bed_filename=self.filename_root + ".bed",
        )

        new = pd.read_csv(
            self.filename_root + ".bed",
            sep="\t",
            header=None,
            names=["contig", "start", "end", "fraction", "read_count"],
        )
        old = pd.read_csv(
            "tests/HG00157-HLA-sorted-test.bed",
            sep="\t",
            header=None,
            names=["contig", "start", "end", "fraction", "read_count"],
        )

        self.assertTrue(
            new[["contig", "start", "end"]].equals(old[["contig", "start", "end"]])
        )


if __name__ == "__main__":
    unittest.main()
