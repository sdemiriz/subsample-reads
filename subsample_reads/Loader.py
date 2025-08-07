import os
import logging
from pathlib import Path
from typing import Optional, Generator
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm

from subsample_reads.Intervals import Intervals
from subsample_reads.FileHandler import FileHandler
from subsample_reads import (
    BATCH_SIZE,
    WRITE_BATCH_SIZE,
    CACHE_SIZE_LIMIT,
    CHUNK_SIZE,
    MAX_WORKERS,
)

logger = logging.getLogger(__name__)


class Loader(FileHandler):

    def __init__(
        self, bam_path: str, template: Optional[pysam.AlignmentFile] = None
    ) -> None:
        """
        Constructor
        """
        logger.info("Loader - Initialize")

        if template:
            super().check_file_exists(path=template.filename)
        self.template = template

        self.bam_path = bam_path
        self.bam = self.load_bam(bam_path=self.bam_path)

        # Add caching for performance
        self._contig_cache = {}
        self._overlap_cache = {}

        logger.info("Loader - Initialized")

    def load_bam(self, bam_path: str) -> pysam.AlignmentFile:
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """
        if self.template:
            logger.info(
                f"Loader - Template file supplied: load BAM {self.bam_path} (write)"
            )
            super().check_file_exists(path=bam_path)
            bam = pysam.AlignmentFile(self.bam_path, mode="wb", template=self.template)

        else:
            logger.info(f"Loader - No template file: load BAM {self.bam_path} (read)")
            super().check_file_exists(path=bam_path)
            bam = pysam.AlignmentFile(self.bam_path, mode="rb")
            bam = self.index_if_needed(bam=bam, bam_path=bam_path)

        return bam

    def index_if_needed(
        self, bam: pysam.AlignmentFile, bam_path: str
    ) -> pysam.AlignmentFile:
        """
        Index the BAM file if it is not already indexed
        """
        try:
            bam.check_index()
        except ValueError as e:
            logger.warning(f"Loader - Indexing BAM file {bam_path}")
            pysam.samtools.index(bam_path)

        return bam

    def run_sampling(
        self,
        bed_dir: str,
        bed_file: str,
        main_seed: int,
        out_bam: str,
        hlala_mode: bool = False,
        hlala_dir: str = "HLA-LA/",
    ) -> None:
        """Sampling method for both regular and HLA*LA modes."""
        logger.info(f"Loader - Begin {'HLA*LA' if hlala_mode else 'regular'} sampling")
        self.main_seed = int(main_seed)
        self.out_bam = out_bam

        self.intervals = self.get_intervals(bed_dir=bed_dir, bed_file=bed_file)
        self.seeds = self.get_interval_seeds(main_seed=self.main_seed)
        self.buckets = self.setup_buckets()

        region_start, region_end = self.intervals.start, self.intervals.end

        if hlala_mode:
            # HLA*LA-specific setup
            self.setup_mapback(hlala_dir=hlala_dir)

            overhang = 1000
            interval_range = (
                region_start - overhang,
                region_end + overhang,
            )

            prg_reads = self.get_prg_reads(contigs=self.get_prg_contigs())
            for r in tqdm(
                iterable=prg_reads, desc="Processed", unit=" PRG-mapped reads"
            ):

                r_chr6 = self.map_read_to_chr6(read=r)
                read_coords = (
                    r_chr6.reference_start,
                    r_chr6.reference_end,
                )

                if self.overlap(
                    read_coords=read_coords,
                    int_coords=interval_range,
                ):
                    self.add_read_to_bucket(read=r_chr6)

        else:
            mapped_reads = self.get_mapped_reads(start=region_start, end=region_end)
            # Use batch processing for better performance
        self.process_reads_in_batches(mapped_reads, batch_size=BATCH_SIZE)

        self.sample_reads_from_buckets()
        self.write_reads()

    def setup_mapback(self, hlala_dir: str) -> None:
        """
        Set up HLA*LA-specific variables
        """

        self.prg_coords_cache = {}

        # GENCODE v48 (non-lncRNA ones selected)
        self.gene_maps = {
            "A": ("chr6", 29941260, 29949572),
            "B": ("chr6", 31353872, 31367067),
            "C": ("chr6", 31268749, 31272130),
            "DMA": ("chr6", 32948613, 32969094),
            "DMB": ("chr6", 32934629, 32941028),
            "DOA": ("chr6", 33004182, 33009591),
            "DPA1": ("chr6", 33064569, 33080775),
            "DPB1": ("chr6", 33075936, 33089696),
            "DQA1": ("chr6", 32628179, 32647062),
            "DQB1": ("chr6", 32659467, 32668383),
            "DRA": ("chr6", 32439878, 32445046),
            "DRB1": ("chr6", 32577902, 32589848),
            "DRB3": ("chr6_GL000250v2_alt", 3824514, 3837642),
            "DRB4": ("chr6_GL000253v2_alt", 3840435, 3855431),
            "E": ("chr6", 30489509, 30494194),
            "F": ("chr6", 29722775, 29738528),
            "G": ("chr6", 29826967, 29831125),
            "H": ("chr6", 29887752, 29890482),
            "K": ("chr6", 29926459, 29929232),
            "L": ("chr6", 30259625, 30261703),
            "MICA": ("chr6", 31399784, 31415315),
            "MICB": ("chr6", 31494881, 31511124),
            "P": ("chr6", 29800415, 29802425),
            "TAP1": ("chr6", 32845209, 32853816),
            "TAP2": ("chr6", 32821833, 32838739),
            "V": ("chr6", 29792234, 29793136),
        }

        # https://github.com/DiltheyLab/ContigAnalysisScripts/blob/master/fasta2svg.py
        self.contig_names = {
            "apd": "chr6_GL000250v2_alt",
            "cox": "chr6_GL000251v2_alt",
            "dbb": "chr6_GL000252v2_alt",
            "mann": "chr6_GL000253v2_alt",
            "mcf": "chr6_GL000254v2_alt",
            "qbl": "chr6_GL000255v2_alt",
            "ssto": "chr6_GL000256v2_alt",
            "chr6": "chr6",
        }

        # From UCSC Browser
        self.alt_contig_maps = {
            "chr6_GL000250v2_alt": ("chr6", 28734408, 33367716),
            "chr6_GL000251v2_alt": ("chr6", 28510120, 33383765),
            "chr6_GL000252v2_alt": ("chr6", 28734408, 33361299),
            "chr6_GL000253v2_alt": ("chr6", 28734408, 33258200),
            "chr6_GL000254v2_alt": ("chr6", 28734408, 33391865),
            "chr6_GL000255v2_alt": ("chr6", 28734408, 33411973),
            "chr6_GL000256v2_alt": ("chr6", 28691466, 33480577),
            "chr6": ("chr6", 1, 33480577),
        }

        sequences_path = (
            Path(hlala_dir) / "graphs" / "PRG_MHC_GRCh38_withIMGT" / "sequences.txt"
        )
        super().check_file_exists(path=str(sequences_path))
        self.sequence_txt = pd.read_csv(
            filepath_or_buffer=sequences_path,
            sep="\t",
            usecols=["Name", "FASTAID"],  # type: ignore
        )

    def get_intervals(self, bed_dir: str, bed_file: str) -> Intervals:
        """
        Set up Interval instances based on BED-provided coordinates
        """
        logger.info("Loader - Ingest Intervals from BED files")
        return Intervals(bed_dir=bed_dir, bed_file=bed_file)

    def get_interval_seeds(self, main_seed: int) -> np.ndarray:
        """
        Generate a seed per interval provided
        """
        logger.info(f"Loader - Generate random seeds")
        np.random.seed(seed=main_seed)
        return np.random.randint(low=0, high=1_000_000, size=len(self.intervals))

    def setup_buckets(self) -> list[list[pysam.AlignedSegment]]:
        """
        Get an empty read bucket to sort reads from per interval provided
        """
        logger.info("Loader - Set up an empty bucket per interval")
        buckets = [[] for i in range(len(self.intervals))]

        # Create a mapping from interval objects to bucket indices for O(1) lookup
        self.interval_to_bucket = {}
        for i, interval in enumerate(self.intervals.tree):
            self.interval_to_bucket[interval] = i

        return buckets

    def get_mapped_reads(
        self, start: int, end: int
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """
        Yield all mapped reads within limits of BED file
        """
        logger.info("Loader - Fetch mapped reads from supplied region")
        contig = self.normalize_contig(contig=self.intervals.contig)
        for r in self.bam.fetch(
            contig=contig,
            start=start,
            end=end,
        ):
            if r.is_mapped:
                yield r

    def overlap(
        self, read_coords: tuple[int, int], int_coords: tuple[int, int]
    ) -> bool:
        """
        Determine whether the read coordinates overlap the interval coordinates (start-end)
        Read cannot hang over the start of the interval
        """
        # Create cache key for overlap calculations
        cache_key = (read_coords, int_coords)
        if cache_key in self._overlap_cache:
            return self._overlap_cache[cache_key]

        result = max(read_coords[0], int_coords[0]) < min(read_coords[1], int_coords[1])
        self._overlap_cache[cache_key] = result
        return result

    def normalize_contig(self, contig: str) -> str:
        """
        Handle both chrN and N contig names (other formats not supported)
        """
        # Use cache for frequently accessed contigs
        if contig in self._contig_cache:
            return self._contig_cache[contig]

        logger.info(f"Loader - Normalize contig name {contig}")

        if contig in self.bam.references:
            self._contig_cache[contig] = contig
            return contig

        if contig.startswith("chr"):
            normalized_contig = contig[3:]
        else:
            normalized_contig = "chr" + contig

        if normalized_contig not in self.bam.references:
            raise ValueError(f"Cannot auto-detect contig name")

        self._contig_cache[contig] = normalized_contig
        return normalized_contig

    def write_reads(self) -> None:
        """
        Open the output BAM file and write all kept reads
        Sort and index for file for random access in following steps
        """
        logger.info(f"Loader - Write reads to file")

        out_bam = Loader(bam_path=self.out_bam, template=self.bam)

        # Write reads in batches for better performance
        for i in range(0, len(self.reads), WRITE_BATCH_SIZE):
            batch = self.reads[i : i + WRITE_BATCH_SIZE]
            for r in batch:
                out_bam.bam.write(read=r)

        out_bam.close()
        self.sort_and_index()

    def sort_and_index(self) -> None:
        """
        Sort, then index file
        """
        logger.info(f"Loader - Sort and index output BAM file")

        temp_file = "temp.bam"
        pysam.sort(self.out_bam, "-o", temp_file)
        os.rename(src=temp_file, dst=self.out_bam)
        pysam.index(self.out_bam)

    def sample_reads_from_buckets(self):
        """
        Sort reads that overlap with BED intervals into buckets
        """
        logger.info("Loader - Sort reads into buckets")

        for i, bucket in enumerate(self.buckets):
            logger.info(f"Loader - {len(bucket)} reads in bucket {i}")

        # After all reads have been sorted into buckets
        self.reads = []
        for bucket, interval, seed in zip(
            self.buckets, self.intervals.tree, self.seeds
        ):
            logger.info(
                f"Loader - Sampling reads in interval [{interval.begin}-{interval.end}]"
            )

            # Count reads that overhang from previous intervals
            overhang_read_count = sum(
                1
                for prev_read in self.reads
                if self.overlap(
                    read_coords=(prev_read.reference_start, prev_read.reference_end),
                    int_coords=(interval.begin, interval.end),
                )
            )

            # Calculate actual amount of reads to sample
            count = int(interval.data) - overhang_read_count

            logger.info(
                f"Loader - {overhang_read_count} reads overlap from previous interval"
            )
            logger.info(f"Loader - {count} reads need to be sampled")

            # Do the sampling
            try:
                np.random.seed(seed=seed)
                self.reads.extend(np.random.choice(a=bucket, size=count, replace=False))
            except ValueError as e:
                logger.warning(f"Loader - No reads in bucket for interval:\n{e}")

    def map_read_to_chr6(self, read: pysam.AlignedSegment) -> pysam.AlignedSegment:
        """
        Map provided PRG-mapped read back to chr6 and corresponding coordinates
        """
        # Convert to dictionary to modify
        r_dict = read.to_dict()

        # Get initial values
        ref_name, ref_pos = r_dict["ref_name"], r_dict["ref_pos"]
        next_ref_name, next_ref_pos = r_dict["next_ref_name"], r_dict["next_ref_pos"]

        # Convert read contig to chr6 and its start coordinate
        r_dict["ref_name"] = "chr6"
        r_dict["ref_pos"] = self.convert_to_chr6(prg_contig=ref_name, prg_coord=ref_pos)

        # Account for "=" as contig name
        if next_ref_name != "=":
            r_dict["next_ref_name"] = "chr6"

        # Switch next read start coordinate to chr6
        r_dict["next_ref_pos"] = self.convert_to_chr6(
            prg_contig=ref_name, prg_coord=next_ref_pos
        )

        return pysam.AlignedSegment.from_dict(
            sam_dict=r_dict, header=self.add_chr6_to_header()
        )

    def get_prg_reads(
        self, contigs: list[str]
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """
        Yield reads mapped to PRG contigs
        """
        for contig in contigs:
            for r in self.bam.fetch(contig=contig):
                if r.is_mapped:
                    yield r

    def get_prg_contigs(self) -> list[str]:
        """
        Get contigs with names that begin with PRG
        """
        return [contig for contig in self.bam.references if contig.startswith("PRG")]

    def convert_to_chr6(self, prg_contig: str, prg_coord: int) -> str:
        """
        Calculate chr6-based start coordinate of a PRG-mapped read
        """
        chr6_coords = self.get_prg_coords(prg_contig=prg_contig)
        return str((chr6_coords[0]) + int(prg_coord))

    def get_prg_coords(self, prg_contig: str) -> tuple[int, int]:
        """
        Return and cache (start, end) coordinates of PRG contig
        """
        if prg_contig not in self.prg_coords_cache:

            sequence_id = self.sequence_txt[self.sequence_txt["FASTAID"] == prg_contig][
                "Name"
            ].values[0]

            gene_name = sequence_id.split("*")[0]

            # If PRG corresponds to HLA allele (with * notation in name)
            if gene_name in self.gene_maps:
                chrom, start, end = self.gene_maps[gene_name]
            # If PRG corresponds to alt contig
            elif sequence_id in self.contig_names:
                alt_contig = self.contig_names[sequence_id]
                chrom, start, end = self.alt_contig_maps[alt_contig]

            self.prg_coords_cache[prg_contig] = (int(start), int(end))

        return self.prg_coords_cache[prg_contig]

    def add_chr6_to_header(self):
        """
        Add chr6 entry to the SQ section of the input BAM file for read back-mapping
        """
        # Get dictionary representation of BAM header and add chr6 entry
        header_sq = self.bam.header.to_dict()["SQ"]
        header_sq.append({"SN": "chr6", "LN": 170805979})  # hg38

        # Replace the original header from BAM with modified one
        header = self.bam.header.to_dict()
        header["SQ"] = header_sq

        # Return de-dictionarified header
        return pysam.AlignmentHeader.from_dict(header_dict=header)

    def add_read_to_bucket(self, read: pysam.AlignedSegment):
        """
        Use IntervalTree to efficiently find overlapping intervals for read assignment
        """
        read_start, read_end = read.reference_start, read.reference_end

        # Use IntervalTree to find overlapping intervals efficiently
        overlapping_intervals = self.intervals.tree.overlap(read_start, read_end)

        if overlapping_intervals:
            # Filter intervals to match your original overlap behavior
            # Your implementation treats adjacent intervals as non-overlapping
            candidate_buckets = []
            for interval in overlapping_intervals:
                # Apply your original overlap logic to ensure consistency
                if self.overlap(
                    read_coords=(read_start, read_end),
                    int_coords=(interval.begin, interval.end),
                ):
                    candidate_buckets.append(self.interval_to_bucket[interval])

            # Randomly select one bucket to deposit the read
            if candidate_buckets:
                np.random.seed(seed=self.main_seed)
                b = np.random.choice(a=candidate_buckets)
                self.buckets[b].append(read)

    def add_read_to_bucket_batch(self, reads: list[pysam.AlignedSegment]) -> None:
        """
        Process multiple reads at once to reduce overhead
        """
        # Pre-calculate read coordinates for all reads
        read_coords = [(r.reference_start, r.reference_end) for r in reads]

        # Batch interval tree lookups
        for read, (read_start, read_end) in zip(reads, read_coords):
            overlapping_intervals = self.intervals.tree.overlap(read_start, read_end)

            candidate_buckets = []
            for interval in overlapping_intervals:
                if self.overlap(
                    read_coords=(read_start, read_end),
                    int_coords=(interval.begin, interval.end),
                ):
                    candidate_buckets.append(self.interval_to_bucket[interval])

            if candidate_buckets:
                np.random.seed(seed=self.main_seed)
                b = np.random.choice(a=candidate_buckets)
                self.buckets[b].append(read)

    def process_reads_in_batches(self, reads_generator, batch_size: int = 1000) -> None:
        """
        Process reads in batches to improve performance
        """
        batch = []
        for read in reads_generator:
            batch.append(read)
            if len(batch) >= batch_size:
                self.add_read_to_bucket_batch(batch)
                batch = []

        # Process remaining reads
        if batch:
            self.add_read_to_bucket_batch(batch)

    def process_region_parallel(
        self,
        start: int,
        end: int,
        chunk_size: int = CHUNK_SIZE,
        max_workers: int = None,
    ) -> list[pysam.AlignedSegment]:
        """
        Process large regions in parallel chunks
        """
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), MAX_WORKERS)

        chunks = []
        for chunk_start in range(start, end, chunk_size):
            chunk_end = min(chunk_start + chunk_size, end)
            chunks.append((chunk_start, chunk_end))

        all_reads = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(self._process_chunk, chunk_start, chunk_end)
                for chunk_start, chunk_end in chunks
            ]

            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Processing chunks"
            ):
                chunk_reads = future.result()
                all_reads.extend(chunk_reads)

        return all_reads

    def _process_chunk(self, start: int, end: int) -> list[pysam.AlignedSegment]:
        """
        Process a single chunk of the genome
        """
        reads = []
        contig = self.normalize_contig(contig=self.intervals.contig)
        for read in self.bam.fetch(contig=contig, start=start, end=end):
            if read.is_mapped:
                reads.append(read)
        return reads

    def fetch(
        self, names: Optional[list[str]] = None
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """
        Yield all mapped reads within limits of BED file
        """
        logger.info("Loader - Fetch mapped reads from supplied region")
        if names is not None:
            for read in self.bam.fetch():
                if read.query_name in names:
                    yield read
        else:
            yield from self.bam.fetch()

    def get_reference_name(self, reference_id: int) -> str:
        """
        Get reference name from reference ID
        """
        return self.bam.get_reference_name(reference_id)

    def clear_caches(self) -> None:
        """
        Clear caches to prevent memory bloat
        """
        if len(self._contig_cache) > 1000:
            self._contig_cache.clear()
        if len(self._overlap_cache) > CACHE_SIZE_LIMIT:
            self._overlap_cache.clear()

    def close(self) -> None:
        """
        Close BAM file using pysam's internal method and clear caches
        """
        logger.info(f"Loader - Close BAM file")
        self.clear_caches()
        self.bam.close()
