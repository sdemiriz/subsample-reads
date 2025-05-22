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
        map_bam: str,
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

        self.in_bam = Loader(file=in_bam)
        self.map_bam = Loader(file=map_bam)
        self.out_bam = Loader(file=out_bam)
        self.bams = {"in": self.in_bam, "map": self.map_bam, "out": self.out_bam}
        self.colormap = {"in": "#648FFF", "map": "#DC267F", "out": "#FFB000"}

        self.out_plt = out_plt
        self.plot()

        info("Plotter - Complete Plotter")

    def get_pileups(self) -> list:
        """
        Pileup BAMs for the defined region
        """
        info("Plotter - Pileup BAMs")

        pileups = [
            bam.bam.pileup(
                contig=bam.normalize_contig(self.intervals.contig),
                start=self.intervals.start,
                end=self.intervals.end,
            )
            for bam in self.bams.values()
        ]

        info("Plotter - Complete pileup BAMs")
        return pileups

    def get_counts(self, start, end) -> list:
        """
        Pileup BAMs for the defined region
        """
        info("Plotter - Pileup BAMs")

        counts = [
            bam.bam.count(
                contig=bam.normalize_contig(self.intervals.contig),
                start=start,
                end=end,
            )
            for bam in self.bams.values()
        ]

        info("Plotter - Complete pileup BAMs")
        return list(counts)

    def plot(self) -> None:
        """
        Plot provided BAM file pileups
        """
        info("Plotter - Begin plotting")
        fig, ax_line, ax_bar = self.plot_setup()

        self.add_lineplot(ax=ax_line)
        self.add_boundaries(ax=ax_line)

        self.add_fractions(ax=ax_line)
        self.add_barplot(ax=ax_bar)

        self.add_annotations(ax_line=ax_line, ax_bar=ax_bar)

        info("Plotter - Save plot")
        plt.savefig(self.out_plt, dpi=3000)

    @staticmethod
    def plot_setup():
        fig, ax_line = plt.subplots(layout="constrained")
        ax_bar = ax_line.twinx()
        return fig, ax_line, ax_bar

    def add_lineplot(self, ax):
        info("Plotter - Plot pileups")

        for p, l, c in zip(
            self.get_pileups(), self.bams.keys(), self.colormap.values()
        ):
            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )
            ax.plot(
                pileup["coord"], pileup["depth"], label=l, alpha=0.8, color=c, zorder=1
            )

    def add_boundaries(self, ax):
        info("Plotter - Draw boundaries")

        for b in self.boundaries:
            ax.axvline(
                x=b, linestyle="--", linewidth=1.5, color="gray", alpha=0.5, zorder=1
            )

    def add_annotations(self, ax_line, ax_bar):
        info("Plotter - Add axes, title, legend")

        ax_line.ticklabel_format(useOffset=False, style="plain")
        ax_line.set_title(
            f"Coverage across {self.intervals.contig}:{self.intervals.start}-{self.intervals.end}"
        )
        ax_line.set_xlabel("Chromosomal coordinate")
        ax_line.set_ylabel("Depth of coverage")
        ax_bar.set_ylabel("Read count")
        ax_line.legend(loc="upper right")

    def add_fractions(self, ax):
        info("Plotter - Add fraction annotations")
        for row in self.intervals.bed.iterrows():
            ax.text(
                x=(row[1]["start"] + row[1]["end"]) / 2,
                y=-10,
                s=str(row[1]["fraction"])[:5],
                ha="center",
                color="black",
                size="x-small",
                zorder=2,
            )

    def add_barplot(self, ax):
        info("Plotter - Add barplot")
        for row in self.intervals.bed.iterrows():
            counts = self.get_counts(row[1]["start"], row[1]["end"])

            for i, c in zip(range(3), self.colormap.values()):
                ax.bar(
                    x=(row[1]["start"] + row[1]["end"]) / 2
                    - (2 - (i + 1)) * (row[1]["end"] - row[1]["start"]) / 3,
                    height=counts[i],
                    width=250,
                    color=c,
                    alpha=0.5,
                    zorder=0,
                )
