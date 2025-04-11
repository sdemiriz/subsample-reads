import pandas as pd
import matplotlib.pyplot as plt
import pysam
from lib.BAMloader import BAMloader


class BAMplotter:

    def __init__(
        self,
        bam_filenames: list[str],
        contig: int,
        start: int,
        end: int,
        out: str = "bamloader-plot.png",
    ):

        bams = [BAMloader(bam) for bam in bam_filenames]
        for bam in bams:
            if not bam.bam.has_index():
                bam.index()

        contig = str(contig)

        pileups = []
        for bam in bams:
            pileups.append(
                [
                    (p.pos, p.n)
                    for p in bam.bam.pileup(contig=contig, start=start, end=end)
                ]
            )

        fig, ax = plt.subplots(layout="constrained")

        for pileup in pileups:
            ax.plot(
                [p[0] for p in pileup],
                [p[1] for p in pileup],
                label=f"{bam.file.split('/')[-1]}",
            )

        ax.grid(visible=True, linestyle="--", linewidth=1)
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across chr{contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        plt.savefig(out)
