from intervaltree import Interval, IntervalTree
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

    def run_subsampling(
        self, contig: str, tree: IntervalTree, initial_seed: int, out_bam
    ) -> None:
        """
        Trigger the subsampling procedure of the class
        """
        info("Start subsampling procedure")

        seeds = self.get_sampling_seeds(initial_seed=initial_seed, count=len(tree))
        contig = self.handle_contig_name(contig=contig)

        for seed, interval in zip(seeds, tree):
            self.subsample_interval(
                contig=contig,
                interval=interval,
                seed=seed,
                out_bam=out_bam,
            )

        info("Complete subsampling procedure")

        info("Close BAM file IO")
        self.bam.close()
        out_bam.bam.close()

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

    def subsample_interval(
        self,
        contig: int,
        interval: Interval,
        seed: int,
        out_bam,
    ) -> None:
        """
        Subsample reads inside the specified interval region based on provided fraction
        """
        info(f"Subsample interval {interval.begin}-{interval.end}")

        keep_reads, drop_reads = self.sample(
            contig=contig,
            interval=interval,
            seed=seed,
        )

        # Save query names of dropped reads
        info("\tSave dropped read query names")
        self.drop_cache.update([dropped_read.query_name for dropped_read in drop_reads])

        # Write kept reads to out_bam
        info("\tWrite kept read pairs to output BAM")

        for read in keep_reads:
            out_bam.bam.write(read=read)

        info("Complete subsample interval")

    def sample(self, contig: str, interval: Interval, seed: int) -> tuple:
        """
        Get reads to drop and reads to keep in provided interval based on seed
        """
        info(f"\tSample reads from interval")

        info("\tSet seed {seed}")
        np.random.seed(seed)

        reads_count = sum(
            1
            for r in self.non_dropped_reads(
                contig=contig,
                start=interval.begin,
                end=interval.end,
            )
        )
        info(f"\tFetch {reads_count} reads from interval")

        drop_count = self.get_drop_count(
            drop_fraction=interval.data,
            total_reads_count=reads_count,
        )

        # Get indices of reads that need to be dropped
        info("\tSample indices for reads to drop")
        drop_indices = np.random.choice(
            a=np.arange(reads_count),
            size=drop_count,
            replace=False,
        )

        # Reads to keep
        info(f"\tGet reads to keep")
        keep = [
            read
            for i, read in enumerate(
                self.non_dropped_reads(
                    contig=contig,
                    start=interval.begin,
                    end=interval.end,
                )
            )
            if i not in drop_indices
        ]
        info(f"\tKeep {len(keep)} reads")

        # Reads to drop
        info(f"\tGet reads to drop")
        drop = [
            read
            for i, read in enumerate(
                self.non_dropped_reads(
                    contig=contig,
                    start=interval.begin,
                    end=interval.end,
                )
            )
            if i in drop_indices
        ]
        info(f"\tDrop {len(drop)} reads")
        info(
            f"\tKeep {len(keep)} reads + Drop {len(drop)} reads = Total {len(keep)+len(drop)} reads"
        )

        assert (
            len(keep) + len(drop) == reads_count
        ), f"Kept {len(keep)} reads + Dropped {len(drop)} reads != Total {reads_count} reads"

        info(f"\tComplete sample reads from interval")
        return keep, drop

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

    def non_dropped_reads(self, contig: str, start: int, end: int):
        """
        Yield reads from specified interval if not dropped before
        """
        info(f"\tGet non-dropped reads")

        for read in self.bam.fetch(contig=str(contig), start=start, end=end):
            if read.query_name not in self.drop_cache:
                yield read

        info(f"\tComplete get non-dropped reads")

    def close(self) -> None:
        self.bam.close()
