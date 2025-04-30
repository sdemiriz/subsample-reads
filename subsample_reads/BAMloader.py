from subsample_reads.Intervals import Intervals
from logging import info, warning
import pysam, math
import numpy as np


class BAMloader:

    def __init__(self, file: str, template: pysam.AlignmentFile = None) -> None:
        """
        Constructor:
        If file is initialized with a template and optionally, a template file
        """
        info("Initialize BAMloader")

        self.file = file
        self.template = template

        self.bam = self.load_bam()
        self.drop_cache = set()

        info("Complete BAMloader")

    def load_bam(self) -> pysam.AlignmentFile:
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """
        info("Load BAM file")

        if self.template:
            info(f"Template file supplied: load BAM in write mode")
            bam = pysam.AlignmentFile(self.file, mode="wb", template=self.template)

        else:
            info(f"Template file not supplied: load BAM in read mode")
            bam = pysam.AlignmentFile(self.file, mode="rb")

        info("Complete load BAM file")
        return bam

    def run_sampling(self, intervals: str, initial_seed: int, out_bam) -> None:
        """
        Trigger the subsampling procedure of the class
        """
        info("Start sampling procedure")

        bed = Intervals(file=intervals)
        contig, start, end = bed.get_limits()
        contig = self.handle_contig_name(contig=contig)

        seeds = self.get_sampling_seeds(initial_seed=initial_seed, count=len(bed.tree))
        all_reads = self.bam.fetch(contig=contig, start=start, end=end)
        buckets = len(bed.tree) * [dict()]
        for r in all_reads:
            try:
                low, hi = min(r.get_reference_positions()), max(
                    r.get_reference_positions()
                )
            except ValueError:
                continue
            for interval, bucket in zip(bed.tree, buckets):
                if self.overlap((low, hi), (interval.begin, interval.end)):
                    if r.query_name in bucket:
                        bucket[r.query_name] += 1
                    else:
                        bucket[r.query_name] = 1

        drop_list = []
        keep_list = []
        for seed, b in zip(seeds, buckets):
            np.random.seed(seed)

            for r in drop_list:
                if r in b:
                    if b[r] > 1:
                        b[r] -= 1
                    else:
                        b.pop(r, None)
            drop = np.random.choice(
                a=b.keys(),
                p=[i / sum(b.values()) for i in b.values()],
                size=round((1 - interval.data) * len(b)),
                replace=False,
            )
            keep = [k for k in b.keys() if k not in drop]

            drop_list.extend(drop)
            keep_list.extend(keep)
            assert len(drop) + len(keep) == len(b)

        info("Complete sampling procedure")

        for r in all_reads:
            if r.query_name in keep_list:
                self.bam.write(r)

        info("Close BAM file IO")
        self.bam.close()
        out_bam.bam.close()

    @staticmethod
    def overlap(pair_x: tuple[int, int], pair_y: tuple[int, int]) -> bool:
        """
        Return whether the two provided intervals overlap
        """
        return max(pair_x[0], pair_y[0]) < min(pair_x[1], pair_y[1])

    def handle_contig_name(self, contig: str) -> str:
        """
        Handle contig names based on contigs from BA< file
        """
        info("Start handle contig names")

        if contig not in self.bam.references:
            if contig.startswith("chr"):
                info("Contig starts with 'chr', trying 'N' format")
                contig = contig[3:]
            else:
                info("Contig does not start with 'chr', trying 'chrN' format")
                contig = "chr" + contig

        if contig not in self.bam.references:
            raise ValueError("Contig name could not be automatically fixed")

        info("Complete handle contig names")
        return contig

    @staticmethod
    def get_sampling_seeds(initial_seed: int, count: int) -> list[int]:
        """
        Generate a number of integer seeds from initial_seed
        """
        info("Generate seeds")

        np.random.seed(initial_seed)
        seeds = list(np.random.randint(low=1, high=1_000_000, size=count))

        info("Complete generate seeds")
        return seeds

    @staticmethod
    def get_drop_count(drop_fraction: float, total_reads_count: int) -> int:
        """
        Calculate the integer count of reads to be dropped
        """
        info(f"\tGet dropped reads count")

        drop_count = math.ceil((1 - drop_fraction) * total_reads_count)
        assert drop_count >= 0, f"Drop count cannot be negative ({drop_count= })"

        try:
            info(
                f"\tDrop {drop_count} / {total_reads_count} = {drop_count / total_reads_count} of reads"
            )
        except ZeroDivisionError:
            warning(f"\tZero reads in interval")

        info(f"\tComplete get dropped reads count")
        return drop_count

    def close(self) -> None:
        self.bam.close()
