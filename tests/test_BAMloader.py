import unittest, os
from subsample_reads.BAMloader import BAMloader


class TestBAMLoader(unittest.TestCase):

    test_bam_filename = "HG00157-HLA-sorted.bam"
    bam = BAMloader(file=test_bam_filename)

    def testContigNotInReferences(self):

        with self.assertRaises(Exception):
            self.bam.get_length(contig="ImpossibleConfigName")


if __name__ == "__main__":
    unittest.main()
