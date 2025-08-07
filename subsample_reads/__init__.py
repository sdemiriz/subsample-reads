"""
Subsample reads toolkit for BAM file processing
"""

# Performance configuration constants
BATCH_SIZE = 1000  # Number of reads to process in batches
WRITE_BATCH_SIZE = 10000  # Number of reads to write in batches
CACHE_SIZE_LIMIT = 10000  # Maximum cache size before clearing
CHUNK_SIZE = 1000000  # Size of genome chunks for parallel processing
MAX_WORKERS = 8  # Maximum number of parallel workers

from .Plotter import Plotter
from .Loader import Loader
from .Mapper import Mapper
from .Comparator import Comparator

__all__ = ["Plotter", "Loader", "Mapper", "Comparator"]
