from subsample_reads.Intervals import Intervals
from logging import info, warning
import pandas as pd
import numpy as np
import pysam, os


class Loader:

    def __init__(self, file: str, template: pysam.AlignmentFile = None) -> None:
        """
        Constructor
        """
        info("Loader - Initialize Loader")

        self.file = file
        self.template = template
        self.load_bam()

        info("Loader - Complete initialize Loader")

    def load_bam(self) -> None:
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """
        if self.template:
            info(f"Loader - Template file supplied: load BAM {self.file} in write mode")
            self.bam = pysam.AlignmentFile(self.file, mode="wb", template=self.template)

        else:
            info(f"Loader - No template file: load BAM {self.file} in read mode")
            self.bam = pysam.AlignmentFile(self.file, mode="rb")

    def sample(
        self,
        bed_dir: str,
        bed_file: str,
        main_seed: int,
        out_bam: str,
    ) -> None:
        """
        Sample BAM file according to interval data provided
        """
        info("Loader - Begin sampling")
        self.main_seed = int(main_seed)
        self.out_bam = out_bam

        self.get_intervals(bed_dir=bed_dir, bed_file=bed_file)
        self.get_interval_seeds(main_seed=self.main_seed)
        self.get_empty_buckets()

        mapped_reads = self.get_mapped_reads(
            start=self.intervals.start, end=self.intervals.end
        )

        # For all mapped reads
        for r in mapped_reads:
            self.add_read_to_bucket(read=r)

        # Sort reads and write reads
        self.sample_reads_in_buckets()
        self.write_reads()

    def get_intervals(self, bed_dir: str, bed_file: str) -> None:
        """
        Set up Interval instances based on BED-provided coordinates
        """
        info("Loader - Ingest Intervals from BED files")
        self.intervals = Intervals(bed_dir=bed_dir, bed_file=bed_file)

    def get_interval_seeds(self, main_seed: int) -> None:
        """
        Generate a seed per interval provided
        """
        info(f"Loader - Generate random seeds")

        np.random.seed(seed=main_seed)
        self.seeds = np.random.randint(low=0, high=1_000_000, size=len(self.intervals))

    def get_empty_buckets(self) -> None:
        """
        Get an empty read bucket to sort reads from per interval provided
        """
        info("Loader - Set up an empty bucket per interval")
        self.buckets = [[] for i in range(len(self.intervals))]

    def get_mapped_reads(self, start: int, end: int):
        """
        Yield all mapped reads within limits of BED file
        """
        info("Loader - Fetch mapped reads from supplied region")
        for r in self.bam.fetch(
            contig=self.normalize_contig(self.intervals.contig), start=start, end=end
        ):
            if r.is_mapped:
                yield r

    @staticmethod
    def overlap(read_coords: tuple[int, int], int_coords: tuple[int, int]) -> bool:
        """
        Determine whether the read coordinates overlap the interval coordinates (start-end)
        Read cannot hang over the start of the interval
        """
        return max(read_coords[0], int_coords[0]) < min(read_coords[1], int_coords[1])

    def normalize_contig(self, contig) -> str:
        """
        Handle both chrN and N contig names (other formats not supported)
        """
        info(f"Loader - Normalize contig name {contig}")

        contig = str(contig)
        if contig in self.bam.references:
            return contig

        if contig.startswith("chr"):
            contig = contig[3:]
        else:
            contig = "chr" + contig

        if contig not in self.bam.references:
            warning(
                f"Loader - Contig name could not be parsed automatically ({contig})"
            )
            raise ValueError("Cannot auto-detect contig name")

        return contig

    def write_reads(self) -> None:
        """
        Open the output BAM file and write all kept reads
        Sort and index for file for random access in following steps
        """
        info(f"Loader - Write reads to file")

        out_bam = Loader(file=self.out_bam, template=self.bam)
        for r in self.reads:
            out_bam.bam.write(read=r)

        out_bam.close()
        self.sort_and_index()

    def sort_and_index(self) -> None:
        """
        Sort, then index file
        """
        info(f"Loader - Sort and index output BAM file")

        temp_file = "temp.bam"
        pysam.sort(self.out_bam, "-o", temp_file)
        os.rename(src=temp_file, dst=self.out_bam)
        pysam.index(self.out_bam)

    def hlala(
        self,
        hlala_dir: str,
        bed_dir: str,
        bed_file: str,
        main_seed: int,
        out_bam: str,
    ) -> None:
        """
        Special sampling procedure for HLA*LA tool output, behaves identically to regular sampling
        but "un-maps" PRG-mapped reads to their reference chr6 locations before sampling
        """
        info(f"Loader - Begin HLA*LA sampling")

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

        self.main_seed = int(main_seed)
        self.out_bam = out_bam

        self.prg_coords_cache = {}
        self.sequence_txt = pd.read_csv(
            hlala_dir + "/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
            sep="\t",
            usecols=["Name", "FASTAID"],
        )

        contigs = self.get_prg_contigs()
        self.get_intervals(bed_dir=bed_dir, bed_file=bed_file)
        self.get_interval_seeds(main_seed=self.main_seed)
        self.get_empty_buckets()

        limit_extension = 1000
        self.interval_range = (
            self.intervals.start - limit_extension,
            self.intervals.end + limit_extension,
        )

        prg_contigs = self.get_prg_reads(contigs=contigs)

        skipped_read_count = 0
        info("Loader - Iterate over PRG reads")
        for r in prg_contigs:

            chr6_read = self.map_read_to_chr6(read=r)

            # Don't consider read if it doesn't overlap any intervals
            if not self.overlap(
                read_coords=(chr6_read.reference_start, chr6_read.reference_end),
                int_coords=self.interval_range,
            ):
                skipped_read_count += 1
                if skipped_read_count % 1000 == 0:
                    info(f"Loader - Skipped {skipped_read_count} reads")

                continue

            self.add_read_to_bucket(read=r)

        self.sample_reads_in_buckets()
        self.write_reads()

    def sample_reads_in_buckets(self):
        """
        Sort reads that overlap with BED intervals into buckets
        """
        info("Loader - Sort reads into buckets")

        for i, bucket in enumerate(self.buckets):
            info(f"Loader - {len(bucket)} total reads deposited into bucket {i}")

        # After all reads have been sorted into buckets
        self.reads = []
        for bucket, interval, seed in zip(
            self.buckets, self.intervals.tree, self.seeds
        ):
            info(
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

            info(f"Loader - {overhang_read_count} reads overlap from previous interval")
            info(f"Loader - {count} reads need to be sampled")

            # Do the sampling
            try:
                np.random.seed(seed=seed)
                self.reads.extend(np.random.choice(a=bucket, size=count, replace=False))
            except ValueError as e:
                info(f"Loader - No reads in bucket for interval:\n{e}")

    def map_read_to_chr6(self, read):
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

    def get_prg_reads(self, contigs):
        """
        Yield reads mapped to PRG contigs
        """
        for contig in contigs:
            for r in self.bam.fetch(contig=contig):
                if r.is_mapped:
                    yield r

    def get_prg_contigs(self) -> None:
        """
        Get contigs with names that begin with PRG
        """
        return [contig for contig in self.bam.references if contig.startswith("PRG")]

    def convert_to_chr6(self, prg_contig, prg_coord) -> str:
        """
        Calculate chr6-based start coordinate of a PRG-mapped read
        """
        chr6_coords = self.get_prg_coords(prg_contig=prg_contig)
        return str((chr6_coords[0]) + int(prg_coord))

    def get_prg_coords(self, prg_contig) -> tuple[int, int]:
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
        Compare read coordinates with every bucket to find overlap
        If overlap is found with multiple buckets, randomly place read into one
        """
        # Keep a tally of buckets a read can fall into
        candidate_buckets = []
        prior_has_overlap = False
        for i, interval in enumerate(self.intervals.tree):

            # print(f"{read_coords=}\n{interval.begin=} {interval.end=}\n")

            has_overlap = self.overlap(
                read_coords=(
                    read.reference_start,
                    read.reference_end,
                ),
                int_coords=(interval.begin, interval.end),
            )

            # Reads should overlap a number of sequential intervals
            if has_overlap:
                prior_has_overlap = True
                candidate_buckets.append(i)
            # If no more overlapping intervals in sequence, no need to check further
            elif prior_has_overlap:
                break

        # Randomly select one bucket to deposit the read
        if len(candidate_buckets) > 0:
            np.random.seed(seed=self.main_seed)
            b = np.random.choice(a=candidate_buckets)
            self.buckets[b].append(read)

    def close(self) -> None:
        """
        Close BAM file using pysam's internal method
        """
        info(f"Loader - Close BAM file")
        self.bam.close()
