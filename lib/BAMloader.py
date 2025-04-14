import pysam, math, logging
import numpy as np


class BAMloader:

    def __init__(self, file: str, template: pysam.AlignmentFile = None):
        """
        Constructor:
        If file is initialized with a template and optionally, a template file
        """

        self.file = file
        logging.info(f"Initialize BAMloader using file {self.file}")

        self.load_bam(template=template)
        self.confirm_index(template=template)

        self.drop_memory = set()
        logging.info(f"BAMloader initialized")

    def load_bam(self, template=None):
        """
        Open in write mode if a template has been provided, otherwise open in read mode
        """
        if template:
            self.bam = pysam.AlignmentFile(self.file, mode="wb", template=template)
            logging.info(
                f"Initialize output BAMloader from file {self.file} and template {template}"
            )
        else:
            self.bam = pysam.AlignmentFile(self.file, mode="rb")
            logging.info(f"Initialize input BAMloader using file {self.file}")

    def get_length(self, contig: int):
        """
        Get length of provided contig in number of base pairs
        """
        logging.info(f"Get length for contig {contig}")

        contig = str(contig)
        assert contig in self.bam.references, f"Given {contig=} not in BAM references"

        return self.bam.lengths[self.bam.references.index(contig)]

    def downsample_reads(
        self,
        out_bam,
        contig: int,
        start: int,
        end: int,
        fraction: float,
        seed: int,
    ):
        """
        Subsample reads inside the specified interval region based on provided fraction
        """
        current_interval = f"\tInterval {contig}:{start}-{end}:"
        logging.info(f"{current_interval} Subsample {fraction} of reads, seed {seed}")

        keep_reads, drop_reads = self.sample(
            contig=contig,
            start=start,
            end=end,
            fraction=fraction,
            seed=seed,
            current_interval=current_interval,
        )

        # Save query names of dropped reads
        self.drop_memory.update(
            [dropped_read.query_name for dropped_read in drop_reads]
        )
        logging.info(f"{current_interval} Save dropped read query names")

        # Write kept reads to out_bam
        for keep_read in keep_reads:
            out_bam.bam.write(read=keep_read)
        logging.info(f"{current_interval} Write kept read pairs to output BAM")

    def sample(self, contig, start, end, fraction, seed, current_interval):
        """"""
        np.random.seed(seed)
        reads_count = sum(
            1 for r in self.non_dropped_reads(contig=contig, start=start, end=end)
        )
        logging.info(f"{current_interval} Fetch {reads_count} from interval")

        drop_count = math.ceil((1.0 - fraction) * reads_count)
        assert drop_count >= 0, f"Drop count negative {drop_count}"

        try:
            logging.info(
                f"{current_interval} Drop {drop_count} / {reads_count} = {drop_count / reads_count} of reads"
            )
        except ZeroDivisionError:
            logging.warning(f"{current_interval} Zero reads in interval")

        # Get indices of reads that need to be dropped
        drop_indices = np.random.choice(
            a=np.arange(reads_count),
            size=drop_count,
            replace=False,
        )

        # Reads to keep
        keep = [
            read
            for i, read in enumerate(
                self.non_dropped_reads(contig=contig, start=start, end=end)
            )
            if i not in drop_indices
        ]
        logging.info(f"{current_interval} Keep {len(keep)} reads")

        # Reads to drop
        drop = [
            read
            for i, read in enumerate(
                self.non_dropped_reads(contig=contig, start=start, end=end)
            )
            if i in drop_indices
        ]
        logging.info(f"{current_interval} Drop {len(drop)} reads")

        logging.info(
            f"{current_interval} Keep {len(keep)} reads + Drop {len(drop)} reads = Total {len(keep)+len(drop)} reads"
        )
        assert (
            len(keep) + len(drop) == reads_count
        ), f"Kept {len(keep)} reads + Dropped {len(drop)} reads != Total {reads_count} reads"

        return keep, drop

    def non_dropped_reads(self, contig, start, end):
        """
        Yield reads from specified interval if not dropped before
        """
        for read in self.bam.fetch(contig=str(contig), start=start, end=end):
            if read.query_name not in self.drop_memory:
                yield read

    def num_reads_in_interval(self, contig, start, end):
        """ """
        return sum(1 for r in self.bam.fetch(contig=str(contig), start=start, end=end))

    def confirm_index(self, template):
        """
        Index file if one doesn't exist.
        Do not index file if opening in write mode (if template provided).
        """
        if not self.bam.has_index() and not template:
            logging.warning(f"No index found, indexing BAM {self.file}")
            pysam.index(self.file)

    def close(self):
        self.bam.close()
