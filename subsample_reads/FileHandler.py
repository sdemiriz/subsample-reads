from pathlib import Path
import logging
from typing import Union

logger = logging.getLogger(__name__)


class FileHandler:
    """Base class for file handling operations."""

    @staticmethod
    def check_file_exists(path: Union[str, Path]) -> None:
        """
        Check if a file exists at the given path. Raise FileNotFoundError if not.

        Args:
            path: Path to the file to check.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p.absolute()}")
        if not p.is_file():
            raise FileNotFoundError(f"Path exists but is not a file: {p.absolute()}")
