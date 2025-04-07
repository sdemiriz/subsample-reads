import unittest, os
from lib.BAMloader import BAMloader


class TestSomething(unittest.TestCase):

    test_bam = "tests/HG002.hs37d5.30x.chr6.bam"
    bam = BAMloader(in_file=test_bam)

    def testContigNotInReferences(self):
        with self.assertRaises(Exception):

            self.bam.get_length(contig="NotGreatConfigName")

    def testIndexCreatedWhenNotPresent(self):
        try:
            self.bam.get_reads(contig=6, start=0, end=100_000)
            os.remove(self.test_bam + ".bai")

        except Exception:
            self.fail("Index not created when it should have")


if __name__ == "__main__":
    unittest.main()
