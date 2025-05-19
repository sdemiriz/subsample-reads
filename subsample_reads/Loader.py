from subsample_reads.Intervals import Intervals
from logging import info, warning
from pathlib import Path
import numpy as np
import pysam, os


class Loader:

    def __init__(self, path: str, template: pysam.AlignmentFile = None) -> None:
        """
        Constructor
        """
        info(f"Loader - Initialize Loader for {path}")

        self.path = path
        self.template = template
        self.load_bam()

        info(f"Loader - Complete initialize Loader for {path}")

    def load_bam(self) -> None:
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """
        if self.template:
            info("Loader - Template file supplied: load BAM in write mode")
            self.bam = pysam.AlignmentFile(
                str(self.path), mode="wb", template=self.template
            )

        else:
            info("Loader - Template file not supplied: load BAM in read mode")
            self.bam = pysam.AlignmentFile(str(self.path), mode="rb")

    def sample(
        self,
        intervals: Intervals,
        main_seed: int,
        out_bam: str,
    ) -> None:
        """
        Sample BAM file according to interval data provided
        """
        info(f"Loader - Begin sampling")

        # Get interval info
        self.intervals = intervals

        # Get multiple seeds for per-bucket randomness
        self.get_interval_seeds(main_seed=main_seed)

        # Get empty read buckets for each interval
        self.get_empty_buckets()

        # For all mapped reads
        for r in self.get_mapped_reads(
            start=self.intervals.start, end=self.intervals.end
        ):

            # Keep a tally of buckets a read can fall into
            candidate_buckets = []
            prior_has_overlap = False
            for i, interval in enumerate(self.intervals.tree):

                has_overlap = self.overlap(
                    (r.reference_start, r.reference_end), (interval.begin, interval.end)
                )

                # Reads should overlap a number of sequential intervals
                if has_overlap:
                    prior_has_overlap = True
                    candidate_buckets.append(i)
                # If no more overlapping intervals in sequence, no need to check further
                elif prior_has_overlap:
                    break

            # Randomly select one bucket to deposit the read
            np.random.seed(seed=main_seed)
            b = np.random.choice(a=candidate_buckets)
            self.buckets[b].append(r)

        # After all reads have been sorted into buckets
        self.reads = []
        for bucket, interval, seed in zip(
            self.buckets, self.intervals.tree, self.seeds
        ):
            np.random.seed(seed=seed)
            self.reads.extend(
                np.random.choice(a=bucket, size=int(interval.data), replace=False)
            )

        # Write kept reads
        self.write_reads(path=out_bam)

    def get_interval_seeds(self, main_seed: int) -> None:
        """
        Generate a seed per interval provided
        """
        info(f"Loader - Generate random seeds")

        np.random.seed(seed=main_seed)
        self.seeds = np.random.randint(low=0, high=1_000_000, size=len(self.intervals))

    def get_empty_buckets(self) -> None:
        """ """
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
    def overlap(pair_x: tuple[int, int], pair_y: tuple[int, int]) -> bool:
        """
        Determine whether two provided 1D intervals overlap
        """
        return max(pair_x[0], pair_y[0]) < min(pair_x[1], pair_y[1])

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

    def write_reads(self, path: Path) -> None:
        """
        Open the output BAM file and write all kept reads
        Sort and index for file for random access in following steps
        """
        info(f"Loader - Write reads to file")

        out_bam = Loader(path=path, template=self.bam)

        for r in self.reads:
            out_bam.bam.write(read=r)

        out_bam.close()
        self.sort_and_index(path=path)

    def sort_and_index(self, path: Path) -> None:
        """
        Sort, then index file
        """
        info(f"Loader - Sort and index output BAM file")

        temp_file = "temp.bam"
        pysam.sort(path, "-o", temp_file)
        os.rename(src=temp_file, dst=path)
        pysam.index(path)

    def close(self) -> None:
        """
        Close BAM file using pysam's internal method
        """
        info(f"Loader - Close BAM file")
        self.bam.close()
