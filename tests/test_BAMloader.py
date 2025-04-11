import unittest, os
from lib.BAMloader import BAMloader


class TestBAMLoader(unittest.TestCase):

    test_bam_filename = "tests/HG002.hs37d5.30x.chr6.bam"
    bam = BAMloader(file=test_bam_filename)

    def testContigNotInReferences(self):

        with self.assertRaises(Exception):
            self.bam.get_length(contig="ImpossibleConfigName")

    def testIndexCreatedWhenNotPresent(self):

        self.assertTrue(os.path.isfile(self.test_bam_filename + ".bai"))
        os.remove(self.test_bam_filename + ".bai")


if __name__ == "__main__":
    unittest.main()
