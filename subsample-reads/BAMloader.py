from intervaltree import Interval, IntervalTree
from logging import info, warning
import pysam, math
import numpy as np


class BAMloader:

    def __init__(self, file: str, template: pysam.AlignmentFile = None):
        """
        Constructor:
        If file is initialized with a template and optionally, a template file
        """
        info(f"Initialize BAMloader from {file}, {template=}")

        self.file = file
        self.load_bam(template=template)
        self.confirm_index(template=template)
        self.drop_cache = set()

        info("Complete BAMloader")

    def load_bam(self, template=None):
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """
        info(f"Load BAM file from {self.file}")

        if template:
            info(f"Template file supplied: load BAM in write mode")
            self.bam = pysam.AlignmentFile(self.file, mode="wb", template=template)

        else:
            info(f"Template file not supplied: load BAM in read mode")
            self.bam = pysam.AlignmentFile(self.file, mode="rb")

        info(f"Complete load BAM file from {self.file}")

    def run_subsampling(
        self, contig: str, tree: IntervalTree, initial_seed: int, out_bam: str
    ):
        """
        Triggers the subsampling procedure of the class
        """
        info("Start subsampling procedure")

        seeds = self.get_sampling_seeds(initial_seed=initial_seed, count=len(tree))
        for seed, interval in zip(seeds, tree):
            self.subsample_interval(
                out_bam=out_bam,
                contig=contig,
                interval=interval,
                seed=seed,
            )

        info("Complete subsampling procedure")

    def handle_contig_name(self, contig: str):
        """
        Handle contig names of the formats 'chrN' or 'N'
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
    def get_sampling_seeds(initial_seed: int, count: int):
        """
        Generate a supplied count of integer seeds using initial_seed
        """
        info(f"Generate seeds")

        np.random.seed(initial_seed)
        seeds = list(np.random.randint(low=1, high=1_000_000, size=count))

        info(f"Complete generate seeds")
        return seeds

    def subsample_interval(
        self,
        out_bam,
        contig: int,
        interval: Interval,
        seed: int,
    ):
        """
        Subsample reads inside the specified interval region based on provided fraction
        """
        info(f"Start subsample interval")

        current_interval = f"\tInterval {contig}:{interval.begin}-{interval.end}:"
        contig = self.handle_contig_name(contig)

        keep_reads, drop_reads = self.sample(
            contig=contig,
            interval=interval,
            seed=seed,
            current_interval=current_interval,
        )

        # Save query names of dropped reads
        info(f"{current_interval} Save dropped read query names")
        self.drop_cache.update([dropped_read.query_name for dropped_read in drop_reads])

        # Write kept reads to out_bam
        info(f"{current_interval} Write kept read pairs to output BAM")
        for read in keep_reads:
            out_bam.bam.write(read=read)

        info(f"Complete subsample interval")

    def sample(self, contig: str, interval: Interval, seed: int, current_interval: str):
        """ """
        info(f"{current_interval} Sample reads from interval")

        info(f"{current_interval} Set seed {seed}")
        np.random.seed(seed)

        reads_count = sum(
            1
            for r in self.non_dropped_reads(
                contig=contig,
                start=interval.begin,
                end=interval.end,
                current_interval=current_interval,
            )
        )
        info(f"{current_interval} Fetch {reads_count} from interval")

        drop_count = self.get_drop_count(
            drop_fraction=interval.data,
            total_reads_count=reads_count,
            current_interval=current_interval,
        )

        # Get indices of reads that need to be dropped
        info(f"{current_interval} Sample indices for reads to drop")
        drop_indices = np.random.choice(
            a=np.arange(reads_count),
            size=drop_count,
            replace=False,
        )

        # Reads to keep
        info(f"{current_interval} Get reads to keep")
        keep = [
            read
            for i, read in enumerate(
                self.non_dropped_reads(
                    contig=contig,
                    start=interval.begin,
                    end=interval.end,
                    current_interval=current_interval,
                )
            )
            if i not in drop_indices
        ]
        info(f"{current_interval} Keep {len(keep)} reads")

        # Reads to drop
        info(f"{current_interval} Get reads to drop")
        drop = [
            read
            for i, read in enumerate(
                self.non_dropped_reads(
                    contig=contig,
                    start=interval.begin,
                    end=interval.end,
                    current_interval=current_interval,
                )
            )
            if i in drop_indices
        ]
        info(f"{current_interval} Drop {len(drop)} reads")
        info(
            f"{current_interval} Keep {len(keep)} reads + Drop {len(drop)} reads = Total {len(keep)+len(drop)} reads"
        )

        assert (
            len(keep) + len(drop) == reads_count
        ), f"Kept {len(keep)} reads + Dropped {len(drop)} reads != Total {reads_count} reads"

        info(f"{current_interval} Complete sample reads from interval")
        return keep, drop

    @staticmethod
    def get_drop_count(
        drop_fraction: float, total_reads_count: int, current_interval: str
    ):
        """
        Calculate the integer count of reads to be dropped
        """
        info(f"{current_interval} Get dropped reads count")

        drop_count = math.ceil((1 - drop_fraction) * total_reads_count)
        assert drop_count >= 0, f"Drop count cannot be negative ({drop_count= })"

        try:
            info(
                f"{current_interval} Drop {drop_count} / {total_reads_count} = {drop_count / total_reads_count} of reads"
            )
        except ZeroDivisionError:
            warning(f"{current_interval} Zero reads in interval")

        info(f"{current_interval} Complete get dropped reads count")
        return drop_count

    def non_dropped_reads(
        self, contig: str, start: int, end: int, current_interval: str
    ):
        """
        Yield reads from specified interval if not dropped before
        """
        info(f"{current_interval} Get non-dropped reads")

        for read in self.bam.fetch(contig=str(contig), start=start, end=end):
            if read.query_name not in self.drop_cache:
                yield read

        info(f"{current_interval} Complete get non-dropped reads")

    def confirm_index(self, template):
        """
        Index file if one doesn't exist, do not index file if opening in write mode
        """
        info(f"Confirm index present")

        if not self.bam.has_index() and not template:
            warning(f"No index found, indexing BAM {self.file}")
            pysam.index(self.file)

        info(f"Complete confirm index present")

    def close(self):
        self.bam.close()
