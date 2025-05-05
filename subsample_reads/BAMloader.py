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

    def run_sampling(self, intervals: str, initial_seed: int, out_bam: str) -> None:
        """
        Trigger the subsampling procedure of the class
        """
        info("Loader - Start sampling all intervals")

        self.bed = Intervals(file=intervals)
        seeds = self.get_sampling_seeds(
            initial_seed=initial_seed, count=len(self.bed.tree)
        )
        self.sample(intervals=intervals, seeds=seeds)

        info("Loader - Complete sampling all intervals")

        self.write_kept_reads(out_bam=out_bam)

        info("Loader - Close input BAM file IO")
        self.bam.close()

    def sample(self, intervals: str, seeds: list):
        """ """
        buckets = self.form_buckets()

        info("Loader - Sample reads from input file")

        self.global_drop = []
        self.global_keep = []

        for seed, bucket, interval in zip(seeds, buckets, self.bed.tree):
            np.random.seed(seed)

            for dropped_read in self.global_drop:
                query_name = dropped_read.query_name

                if query_name in self.bucket:
                    if bucket[query_name] > 1:
                        bucket[query_name] -= 1
                    else:
                        bucket.pop(query_name, None)

            drop = np.random.choice(
                a=bucket.keys(),
                p=[i / sum(bucket.values()) for i in bucket.values()],
                size=round((1 - interval.data) * len(bucket)),
                replace=False,
            )
            keep = [k for k in bucket.keys() if k not in drop]

            assert len(drop) + len(keep) == len(bucket)

            self.global_drop.extend(drop)
            self.global_keep.extend(keep)

            info(
                f"Loader - Complete process interval {self.bed.contig}:{interval.begin}-{interval.end}"
            )

        info("Loader - Complete sample reads from input file")

    def write_kept_reads(self, out_bam):
        """ """
        info("Loader - Write reads to output file")

        out_bam = BAMloader(file=out_bam, template=self.bam)
        start, end = self.bed.get_limits()

        for r in self.bam.fetch(
            contig=self.handle_contig_name(contig=self.bed.contig), start=start, end=end
        ):
            if r.query_name in self.global_keep:
                out_bam.bam.write(r)

        info("Loader - Complete write reads to output file")
        out_bam.bam.close()
        self.index_bam(out_bam.file)

    @staticmethod
    def index_bam(filename: str) -> None:
        """ """
        info("Loader - Index output BAM")

        pysam.index(filename)

        info("Loader - Complete index output BAM")

    def form_buckets(self):
        """ """
        info("Loader - Form buckets from provided intervals")

        start, end = self.bed.get_limits()
        buckets = len(self.bed.tree) * [dict()]

        for r in self.bam.fetch(
            contig=self.handle_contig_name(contig=self.bed.contig), start=start, end=end
        ):
            if r.is_mapped:
                pos = r.get_reference_positions()
                low, hi = min(pos), max(pos)

                for interval, bucket in zip(self.bed.tree, buckets):
                    if self.overlap((low, hi), (interval.begin, interval.end)):
                        if r.query_name in bucket:
                            bucket[r.query_name] += 1
                        else:
                            bucket[r.query_name] = 1

        info("Loader - Complete form buckets from provided intervals")
        print(len(buckets))
        return buckets

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
        info("Loader - Start handle contig names")

        if contig not in self.bam.references:
            if contig.startswith("chr"):
                info("Loader - Contig starts with 'chr', trying 'N' format")
                contig = contig[3:]
            else:
                info("Loader - Contig does not start with 'chr', trying 'chrN' format")
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
        info("Loader - Generate seeds")

        np.random.seed(initial_seed)
        seeds = list(np.random.randint(low=1, high=1_000_000, size=count))

        info("Loader - Complete generate seeds")
        return seeds

    @staticmethod
    def get_drop_count(drop_fraction: float, total_reads_count: int) -> int:
        """
        Calculate the integer count of reads to be dropped
        """
        info("Loader - Get dropped reads count")

        drop_count = math.ceil((1 - drop_fraction) * total_reads_count)
        assert drop_count >= 0, f"Drop count cannot be negative ({drop_count= })"

        try:
            info(
                f"Loader - Drop {drop_count} / {total_reads_count} = {drop_count / total_reads_count} of reads"
            )
        except ZeroDivisionError:
            warning(f"Loader - Zero reads in interval")

        info("Loader - Complete get dropped reads count")
        return drop_count

    def close(self) -> None:
        self.bam.close()
