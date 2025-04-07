import pysam
import matplotlib.pyplot as plt
from collections import defaultdict


class BAMloader:

    def __init__(self, in_file: str):

        self.in_file = in_file
        self.load_bam()

        self.read_dict = defaultdict(lambda: [None, None])

    def load_bam(self):
        self.bam = pysam.AlignmentFile(self.in_file, "rb")

    def get_length(self, contig: int):

        contig = str(contig)
        assert contig in self.bam.references, f"Given {contig=} not in BAM references"

        return self.bam.lengths[self.bam.references.index(contig)]

    def get_read_pairs(self, contig: int, start: int, end: int):

        for read in self.bam.fetch(contig=str(contig), start=start, end=end):
            if not read.is_proper_pair:
                continue

            qname = read.query_name
            if qname not in self.read_dict:
                if read.is_read1:
                    self.read_dict[qname][0] = read
                else:
                    self.read_dict[qname][1] = read

            else:
                if read.is_read1:
                    yield read, self.read_dict[qname][1]
                else:
                    yield self.read_dict[qname][0], read
                del self.read_dict[qname]

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
            label=f"{self.in_file.split('/')[-1]}",
        )

        ax.grid(visible=True, linestyle="--", linewidth=1)
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across chr{contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        plt.savefig(out)
