import pysam, math
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict


class BAMloader:

    def __init__(self, file: str, template=None):
        self.file = file
        self.load_bam(template=template)

        if not self.bam.has_index():
            print(f"Warning: No index found, indexing BAM {self.file}")
            pysam.index(self.file)

        self.read_dict = defaultdict(lambda: [None, None])
        self.dropped_read_pairs = list()

    def load_bam(self, template=None):

        if template:
            self.bam = pysam.AlignmentFile(self.file, mode="wb", template=template)
        else:
            self.bam = pysam.AlignmentFile(self.file, mode="rb")

    def get_length(self, contig: int):

        contig = str(contig)
        assert contig in self.bam.references, f"Given {contig=} not in BAM references"

        return self.bam.lengths[self.bam.references.index(contig)]

    def get_read_pairs(self, contig: int, start: int, end: int):

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

    def downsample_reads(
        self,
        out_bam,
        contig: int,
        start: int,
        end: int,
        fraction: float,
    ):
        # Get paired reads within the interval
        all_reads_in_interval = list(
            self.get_read_pairs(contig=str(contig), start=start, end=end)
        )

        # Get paired reads that have not been dropped before
        predropped_reads = [
            read_pair
            for read_pair in all_reads_in_interval
            if read_pair[0].query_name in self.dropped_read_pairs
        ]

        # Calculate how many paired reads need to be dropped
        base_size = math.ceil((1.0 - fraction) * len(all_reads_in_interval))

        # Get indices of reads that need to be dropped
        remove_indices = np.random.choice(
            a=len(all_reads_in_interval),
            size=base_size - len(predropped_reads),
            replace=False,
        )

        # Reads to add to out_bam
        kept_reads = [
            read_pair
            for i, read_pair in enumerate(all_reads_in_interval)
            if i not in remove_indices
        ]

        # Reads to remove
        removed_reads = [
            read_pair
            for i, read_pair in enumerate(all_reads_in_interval)
            if i in remove_indices
        ]

        removed_reads += predropped_reads

        # Add removed reads to memory for future intervals
        for removed_read_pair in removed_reads:
            self.dropped_read_pairs.append(removed_read_pair[0].query_name)

        # Write kept reads to out_bam
        for read_pair in kept_reads:
            out_bam.bam.write(read=read_pair[0])
            out_bam.bam.write(read=read_pair[1])

        print(
            f"{str(contig)} {str(start)} {str(end)}, {len(all_reads_in_interval)=}, {len(removed_reads)=}, {len(self.dropped_read_pairs)=}"
        )

    def index(self):
        pysam.index(self.file)

    def plot_pileup(
        self, contig: int, start: int, end: int, out: str = "bamloader-plot.png"
    ):
        contig = str(contig)

        pileup = [
            (p.pos, p.n) for p in self.bam.pileup(contig=contig, start=start, end=end)
        ]

        fig, ax = plt.subplots(layout="constrained")

        ax.plot(
            [p[0] for p in pileup],
            [p[1] for p in pileup],
            label=f"{self.file.split('/')[-1]}",
        )

        ax.grid(visible=True, linestyle="--", linewidth=1)
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across chr{contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        plt.savefig(out)
