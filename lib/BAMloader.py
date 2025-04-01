import pysam
import matplotlib.pyplot as plt


class BAMloader:

    def __init__(self, in_file: str):

        self.in_file = in_file

        self.load_bam()

    def load_bam(self):
        self.bam = pysam.AlignmentFile(self.in_file, "rb")

    def get_reads(self, contig: int, start: int, end: int):
        return self.bam.fetch(contig=str(contig), start=start, end=end)

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
