from subsample_reads.BAMloader import BAMloader
from subsample_reads.Intervals import Intervals
import matplotlib.pyplot as plt
from logging import info
import pandas as pd


class Plotter:

    def __init__(
        self,
        bam_files: list[str],
        bed_file: str,
        out: str | None,
    ) -> None:
        """
        Constructor for plotting utility
        """
        info(f"Plotter - Initialize Plotter with {bam_files=}, {bed_file=}, {out=}")

        self.bed_file = bed_file
        self.bed = Intervals(file=self.bed_file)
        self.boundaries = set(list(self.bed.bed["start"]) + list(self.bed.bed["end"]))

        self.bam_files = bam_files

        self.out = out
        self.plot()

        info("Plotter - Complete Plotter")

    def get_pileups(self) -> list:
        """
        Pileup BAMs for the defined region
        """
        info("Plotter - Pileup BAMs")

        start, end = self.bed.get_limits()
        bams = [BAMloader(file=bam) for bam in self.bam_files]
        pileups = [
            bam.bam.pileup(
                contig=bam.normalize_contig(self.bed.contig), start=start, end=end
            )
            for bam in bams
        ]

        info("Plotter - Complete pileup BAMs")
        return pileups

    def plot(self) -> None:
        """
        Plot provided BAM file pileups
        """
        info("Plotter - Begin plotting")
        fig, ax = plt.subplots(layout="constrained")

        start, end = self.bed.get_limits()
        pileups = self.get_pileups()

        info("Plotter - Iterate pileups")
        for p, b in zip(pileups, self.bam_files):

            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )

            ax.plot(
                pileup["coord"],
                pileup["depth"],
                label=f"{b.split('.')[-2]}",
                alpha=0.5,
            )
        info("Plotter - Complete iterate pileups")

        for b in self.boundaries:
            ax.axvline(x=b, linestyle="--", linewidth=1.5, color="red", alpha=0.3)

        for row in self.bed.bed.iterrows():
            ax.text(
                x=(row[1]["start"] + row[1]["end"]) / 2,
                y=-10,
                s=str(row[1]["fraction"])[:5],
                ha="center",
                size="x-small",
            )

        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across {self.bed.contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        info("Plotter - Complete plotting")

        info("Plotter - Save plot")
        plt.savefig(self.out)
