import pysam, math, logging
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict


class BAMloader:

    def __init__(self, file: str, template: pysam.AlignmentFile = None):
        """
        Constructor:
        If file is initialized with a template and optionally, a template file
        """

        self.file = file
        logging.info(f"Initialized class using file {self.file}")

        self.load_bam(template=template)
        self.confirm_index(template=template)

        self.read_dict = defaultdict(lambda: [None, None])
        self.dropped_read_pairs = list()

    def load_bam(self, template=None):
        """
        Open in write mode if a template has been provided, otherwise open in read mode
        """

        if template:
            self.bam = pysam.AlignmentFile(self.file, mode="wb", template=template)
            logging.info(
                f"Initialized class using file {self.file} and template {template} in write mode"
            )
        else:
            self.bam = pysam.AlignmentFile(self.file, mode="rb")
            logging.info(f"Initialized class using file {self.file} in read mode")

    def get_length(self, contig: int):
        """
        Get length of provided contig in number of base pairs
        """
        contig = str(contig)
        assert contig in self.bam.references, f"Given {contig=} not in BAM references"

        logging.info(f"Getting length for contig {contig}")

        return self.bam.lengths[self.bam.references.index(contig)]

    def get_read_pairs(self, contig: int, start: int, end: int):
        """
        Maintain a dictionary of read pairs using the shared query_name as the key
        """
        logging.info(f"Interval {contig}:{start}-{end}: Fetching read pairs for ")

        for read in self.bam.fetch(contig=str(contig), start=start, end=end):
            if not read.is_proper_pair:
                continue

            qname = read.query_name
            if qname in self.read_dict:
                if read.is_read1:
                    yield read, self.read_dict[qname][1]
                else:
                    yield self.read_dict[qname][0], read

            else:
                if read.is_read1:
                    self.read_dict[qname][0] = read
                else:
                    self.read_dict[qname][1] = read

        logging.info(
            f"{len(self.read_dict)} read pairs for interval {contig}:{start}-{end}"
        )

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
        logging.info(
            f"Interval {contig}:{start}-{end}: Subsampling {fraction} of reads using seed {seed}"
        )
        np.random.seed(seed)

        # Get paired reads within the interval
        all_reads_in_interval = list(
            self.get_read_pairs(contig=str(contig), start=start, end=end)
        )

        # Get paired reads that have been dropped before
        # predropped_reads = [
        #     read_pair
        #     for read_pair in all_reads_in_interval
        #     if read_pair[0].query_name in self.dropped_read_pairs
        # ]

        # Calculate how many paired reads need to be dropped
        base_size = math.ceil((1.0 - fraction) * len(all_reads_in_interval))

        logging.info(
            f"Interval {contig}:{start}-{end}: Dropping a maximum of {base_size} reads out of {len(all_reads_in_interval)}"
        )
        logging.info(
            f"Interval {contig}:{start}-{end}: Actual ratio of reads dropped: {base_size / len(all_reads_in_interval)}"
        )

        # Get indices of reads that need to be dropped
        remove_indices = np.random.choice(
            a=len(all_reads_in_interval),
            size=base_size,  # - len(predropped_reads),
            replace=False,
        )

        # Reads to add to out_bam
        kept_reads = [
            read_pair
            for i, read_pair in enumerate(all_reads_in_interval)
            if i not in remove_indices
        ]

        logging.info(
            f"Interval {contig}:{start}-{end}: Keeping {len(kept_reads)} read pairs"
        )

        # Reads to remove
        removed_reads = [
            read_pair
            for i, read_pair in enumerate(all_reads_in_interval)
            if i in remove_indices
        ]

        logging.info(
            f"Interval {contig}:{start}-{end}: Dropping {len(removed_reads)} read pairs"
        )

        # removed_reads += predropped_reads

        # Add removed reads to memory for future intervals
        # for removed_read_pair in removed_reads:
        #     self.dropped_read_pairs.append(removed_read_pair[0].query_name)

        logging.info(
            f"Interval {contig}:{start}-{end}: Writing kept read pairs to output BAM"
        )
        # Write kept reads to out_bam
        for read_pair in kept_reads:
            out_bam.bam.write(read=read_pair[0])
            out_bam.bam.write(read=read_pair[1])

        logging.info(f"Interval {contig}:{start}-{end}: Subsampling finished.")

    def confirm_index(self, template):
        """
        Index file if one doesn't exist.
        Do not index file if opening in write mode (if template provided)
        """
        if not self.bam.has_index() and not template:
            logging.WARN(f"No index found, indexing BAM {self.file}")
            pysam.index(self.file)
