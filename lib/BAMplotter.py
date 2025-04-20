import pandas as pd
import matplotlib.pyplot as plt
import pysam
from logging import info
from lib.BAMloader import BAMloader
from lib.Intervals import Intervals


class BAMplotter:

    def __init__(
        self,
        bam_files: list[str],
        bed_file: str,
        out: str,
    ):
        self.bed_file = bed_file
        self.bed = Intervals(self.bed_file)
        contig, start, end = self.bed.limits()

        current_interval = f"\tInterval {contig}:{start}-{end}:"
        info(f"{current_interval} Start plot")

        info(f"{current_interval} Load {len(bam_files)} BAM files")
        self.bam_files = bam_files
        self.bams = [BAMloader(bam) for bam in bam_files]

        self.out = out
        self.plot(out=self.out, current_interval=current_interval)

    def plot(self, out, current_interval):
        """
        Plot provided BAM file pileups
        """
        contig, start, end = self.bed.limits()
        info(f"{current_interval} Pileup BAMs")

        pileups = [
            bam.bam.pileup(contig=contig, start=start, stop=end) for bam in self.bams
        ]

        info(f"{current_interval} Begin plot")
        fig, ax = plt.subplots(layout="constrained")

        for pileup in pileups:
            ax.plot(
                [p[0] for p in pileup],
                [p[1] for p in pileup],
                label=f"{self.bam.file.split('/')[-1]}",
            )

        ax.grid(visible=True, linestyle="--", linewidth=1)
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across chr{contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()
        info(f"{current_interval} End plot")

        info(f"{current_interval} Save plot")
        plt.savefig(out)
