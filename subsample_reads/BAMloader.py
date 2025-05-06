from subsample_reads.Intervals import Intervals
from logging import info, warning
import pysam, os
import numpy as np


class BAMloader:

    def __init__(self, file: str, template: pysam.AlignmentFile = None) -> None:
        """
        Constructor:
        If file is initialized with a template and optionally, a template file
        """
        info("Loader - Initialize BAMloader")

        self.file = file
        self.template = template
        self.bam = self.load_bam()

        info("Loader - Complete BAMloader")

    def load_bam(self) -> pysam.AlignmentFile:
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """
        info("Loader - Load BAM file")

        if self.template:

            info("Loader - Template file supplied: load BAM in write mode")
            bam = pysam.AlignmentFile(self.file, mode="wb", template=self.template)

        else:
            info("Loader - Template file not supplied: load BAM in read mode")
            bam = pysam.AlignmentFile(self.file, mode="rb")

        info("Loader - Complete load BAM file")
        return bam

    def sample(
        self,
        intervals: str,
        initial_seed: int,
        out_bam: str,
    ) -> None:

        info(f"Loader - Begin sampling")

        # Get interval info
        self.bed = Intervals(file=intervals)
        start, end = self.bed.get_limits()
        info(f"Loader - Full region {start}-{end}")

        # Get multiple seeds for per-bucket randomness
        self.seeds = self.get_interval_seeds(initial_seed=initial_seed)

        # Get empty read buckets for each interval
        self.buckets = [[] for i in range(len(self.bed.tree))]

        info(f"Loader - Sort reads into buckets")
        # For all reads
        total_read_count = 0
        for r in self.bam.fetch(
            contig=self.normalize_contig(self.bed.contig), start=start, end=end
        ):

            # Given read is mapped
            if r.is_unmapped:
                continue

            # Keep a tally of kept reads
            total_read_count += 1

            # Keep a tally of buckets a read can fall into
            candidate_buckets = []
            for i, interval in enumerate(self.bed.tree):

                if self.overlap(
                    (r.reference_start, r.reference_end), (interval.begin, interval.end)
                ):
                    candidate_buckets.append(i)

            # Randomly select one bucket to deposit the read
            np.random.seed(seed=initial_seed)
            b = np.random.choice(a=candidate_buckets)
            self.buckets[b].append(r)

        info(f"Loader - Complete sort reads into buckets")

        # After all reads have been sorted into buckets
        self.reads = []
        for bucket, interval, seed in zip(self.buckets, self.bed.tree, self.seeds):
            np.random.seed(seed=seed)

            # Sample each bucket for the pre-calculated ratio of reads
            size = round(interval.data * total_read_count)
            if size > len(bucket):
                size = len(bucket)
            bucket = np.random.choice(a=bucket, size=size, replace=False)
            self.reads.extend(bucket)

        # Write kept reads
        self.write_reads(filename=out_bam)

    def get_interval_seeds(self, initial_seed: int):
        """
        Generate a number of reads equal to intervals received from BED file
        """
        np.random.seed(initial_seed)
        info(f"Loader - Generate random seeds")
        return np.random.randint(low=0, high=1000000, size=len(self.bed.tree))

    @staticmethod
    def overlap(pair_x: tuple[int, int], pair_y: tuple[int, int]) -> bool:
        """
        Confirm whether two provided 1D intervals overlap
        """
        return max(pair_x[0], pair_y[0]) < min(pair_x[1], pair_y[1])

    def normalize_contig(self, contig) -> str:
        """
        Handle both chrN and N contig name descriptions
        (other formats not supported)
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

    def write_reads(self, filename: str) -> None:
        """
        Open the output BAM file and write all kept reads
        Sort and index for file for random access in following steps
        """
        info(f"Loader - Write reads to file")
        out_bam = BAMloader(file=filename, template=self.bam)

        for r in self.reads:
            out_bam.bam.write(read=r)

        out_bam.close()
        self.sort_and_index(filename=filename)

    def sort_and_index(self, filename: str) -> None:
        """
        Sort, then index file
        """
        info(f"Loader - Sort and index output BAM file")
        temp_file = "temp.bam"
        pysam.sort(filename, "-o", temp_file)
        os.rename(src=temp_file, dst=filename)
        pysam.index(filename)

    def close(self) -> None:
        """
        Close BAM file using pysam's internal method
        """
        info(f"Loader - Sort and index output BAM file")
        self.bam.close()
