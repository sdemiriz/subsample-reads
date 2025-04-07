import unittest
from lib.BEDloader import BEDloader


class TestBEDloader(unittest.TestCase):

    def testNonUniqueContigsInBED(self):
        with self.assertRaises(Exception):
            BEDloader(in_file="tests/not-unique.bed")


if __name__ == "__main__":
    unittest.main()
