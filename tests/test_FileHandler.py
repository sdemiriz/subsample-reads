import os
import tempfile
import unittest

from subsample_reads.FileHandler import FileHandler


class TestFileHandler(unittest.TestCase):
    def test_check_file_exists_success(self):
        # Create a temporary file and check it exists
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_path = tmp.name
        try:
            FileHandler.check_file_exists(tmp_path)
        finally:
            os.remove(tmp_path)

    def test_check_file_exists_failure(self):
        # Use a filename that should not exist
        with self.assertRaises(FileNotFoundError):
            FileHandler.check_file_exists("DOES_NOT_EXIST.txt")


if __name__ == "__main__":
    unittest.main()
