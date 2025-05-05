from subsample_reads.Intervals import Intervals
from logging import info, warning
import pysam, math, os
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
        self.drop_cache = set()

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
        intervals,
        initial_seed: int,
        out_bam: str,
    ) -> None:

        self.bed = Intervals(file=intervals)
        start, end = self.bed.get_limits()
        self.seeds = self.get_interval_seeds(initial_seed=initial_seed)

        info(f"Loader - Full region {start}-{end}")

        self.reads = []
        for r in self.bam.fetch(
            contig=self.normalize_contig(self.bed.contig), start=start, end=end
        ):
            self.reads.append(r)

        self.write_reads(filename=out_bam)

    def get_interval_seeds(self, initial_seed: int):
        return np.random.randint(low=0, high=1000000, size=len(self.bed.tree))

    @staticmethod
    def overlap(pair_x: tuple[int, int], pair_y: tuple[int, int]):
        return max(pair_x[0], pair_y[0]) < min(pair_x[1], pair_y[1])

    def normalize_contig(self, contig):
        contig = str(contig)

        if contig in self.bam.references:
            return contig

        if contig.startswith("chr"):
            contig = contig[3:]
        else:
            contig = "chr" + contig

        if contig not in self.bam.references:
            raise ValueError("Cannot auto-detect contig name")

        return contig

    def write_reads(self, filename: str) -> None:

        out_bam = BAMloader(file=filename, template=self.bam)

        for r in self.reads:
            out_bam.bam.write(read=r)

        out_bam.close()
        self.sort_and_index(filename=filename)

    def sort_and_index(self, filename: str):
        """ """
        temp_file = "temp.bam"
        pysam.sort(filename, "-o", temp_file)
        pysam.sort(temp_file, "-o", filename)
        os.remove(temp_file)
        pysam.index(filename)

    def close(self) -> None:
        self.bam.close()
