"""
Subsample reads toolkit for BAM file processing

This package provides tools for mapping, sampling, and plotting BAM files
with support for both regular sampling and HLA-LA PRG back-mapping.
"""

# Performance configuration constants
BATCH_SIZE: int = 1000  # Number of reads to process in batches
WRITE_BATCH_SIZE: int = 10000  # Number of reads to write in batches
CACHE_SIZE_LIMIT: int = 10000  # Maximum cache size before clearing
CHUNK_SIZE: int = 1000000  # Size of genome chunks for parallel processing
MAX_WORKERS: int = 8  # Maximum number of parallel workers

from .Plotter import Plotter
from .Loader import Loader
from .Mapper import Mapper
from .Comparator import Comparator

__all__ = ["Plotter", "Loader", "Mapper", "Comparator"]
