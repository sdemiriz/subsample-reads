from subsample_reads.Intervals import Intervals
from subsample_reads.Loader import Loader
import matplotlib.pyplot as plt
from pathlib import Path
from logging import info
import pandas as pd


class Plotter:

    def __init__(
        self,
        in_bam: str,
        out_bam: str,
        bed_dir: str,
        bed_file: str,
        out_plt: str | None,
    ) -> None:
        """
        Constructor for plotting utility
        """
        info(f"Plotter - Initialize Plotter")

        self.intervals = Intervals(bed_dir=bed_dir, bed_file=bed_file)
        self.boundaries = set(
            list(self.intervals.bed["start"]) + list(self.intervals.bed["end"])
        )

        self.in_bam = in_bam
        self.out_bam = out_bam
        self.out_plt = out_plt
        self.plot()

        info("Plotter - Complete Plotter")

    def get_pileups(self) -> list:
        """
        Pileup BAMs for the defined region
        """
        info("Plotter - Pileup BAMs")

        bams = [Loader(file=self.in_bam), Loader(file=self.out_bam)]
        pileups = [
            bam.bam.pileup(
                contig=bam.normalize_contig(self.intervals.contig),
                start=self.intervals.start,
                end=self.intervals.end,
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

        pileups = self.get_pileups()

        info("Plotter - Iterate pileups")
        for p, b in zip(pileups, [self.in_bam, self.out_bam]):

            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )

            ax.plot(
                pileup["coord"],
                pileup["depth"],
                label=Path(b).stem,
                alpha=0.5,
            )
        info("Plotter - Complete iterate pileups")

        for b in self.boundaries:
            ax.axvline(x=b, linestyle="--", linewidth=1.5, color="red", alpha=0.3)

        for row in self.intervals.bed.iterrows():
            ax.text(
                x=(row[1]["start"] + row[1]["end"]) / 2,
                y=-10,
                s=str(row[1]["read_count"])[:5],
                ha="center",
                size="x-small",
            )

        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(
            f"Coverage across {self.intervals.contig}:{self.intervals.start}-{self.intervals.end}"
        )
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        info("Plotter - Complete plotting")

        info("Plotter - Save plot")
        plt.savefig(self.out_plt)
