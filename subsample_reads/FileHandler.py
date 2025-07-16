from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class FileHandler:

    @staticmethod
    def check_file_exists(path: str) -> None:
        """
        Check if a file exists at the given path. Raise FileNotFoundError if not.
        Args:
            path: Path to the file.
        """
        p = Path(path)
        if not p.exists():
            logger.error(f"FileHandler - File not found: {p}")
            raise FileNotFoundError(f"FileHandler - File not found: {p}")
