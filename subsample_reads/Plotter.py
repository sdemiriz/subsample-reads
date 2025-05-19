from subsample_reads.Intervals import Intervals
from subsample_reads.Loader import Loader
import matplotlib.pyplot as plt
from pathlib import Path
from logging import info
import pandas as pd


class Plotter:

    def __init__(
        self,
        in_bam_path: Path,
        map_bam_paths: list[Path],
        out_bam_path: Path,
        intervals: Intervals,
        plt_path: Path,
    ) -> None:
        """
        Constructor for plotting utility
        """
        info(f"Plotter - Initialize Plotter")

        self.intervals = intervals
        self.boundaries = set(
            list(self.intervals.beds[0]["start"]) + list(self.intervals.beds[0]["end"])
        )

        self.in_bam_path = in_bam_path
        self.map_bam_paths = map_bam_paths
        self.out_bam_path = out_bam_path

        self.all_bam_paths = (
            [self.in_bam_path] + self.map_bam_paths + [self.out_bam_path]
        )

        self.plt_path = plt_path
        self.plot()

        info("Plotter - Complete Plotter")

    def plot(self) -> None:
        """
        Plot provided BAM file pileups
        """
        info("Plotter - Begin plotting")
        fig, ax = plt.subplots(layout="constrained")

        pileups = self.get_pileups()

        info("Plotter - Iterate pileups")
        for p, b in zip(pileups, self.all_bam_paths):

            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )

            ax.plot(
                pileup["coord"],
                pileup["depth"],
                label=b.stem,
                alpha=0.5,
            )
        info("Plotter - Complete iterate pileups")

        ax.plot(
            (self.intervals.stats["start"] + self.intervals.stats["end"]) / 2,
            self.intervals.stats["mean"],
            color="b",
        )
        ax.fill_between(
            x=(self.intervals.stats["start"] + self.intervals.stats["end"]) / 2,
            y1=self.intervals.stats["min"],
            y2=self.intervals.stats["max"],
            alpha=0.3,
            color="b",
        )

        for b in self.boundaries:
            ax.axvline(x=b, linestyle="--", linewidth=1.5, color="red", alpha=0.3)

        for row in self.intervals.stats.iterrows():
            ax.text(
                x=(row[1]["start"] + row[1]["end"]) / 2,
                y=-10,
                s=str(row[1]["mean"])[:5],
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
        plt.savefig(self.plt_path)

    def get_pileups(self) -> list:
        """
        Pileup BAMs for the defined region
        """
        info("Plotter - Pileup BAMs")

        bams = [Loader(path=p) for p in self.all_bam_paths]
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
