from subsample_reads.FileHandler import FileHandler
from subsample_reads.Intervals import Intervals
from subsample_reads.Loader import Loader
import matplotlib.pyplot as plt
import pandas as pd
import logging
import pysam

logger = logging.getLogger(__name__)


class Plotter(FileHandler):

    def __init__(
        self,
        in_bam: str,
        map_bam: str,
        sub_bam: str,
        bed_dir: str,
        bed: str,
        out_plt: str,
    ) -> None:
        """
        Initialize Plotter with argument names.

        Args:
            in_bam:   Path to input/original BAM file.
            map_bam:  Path to mapping BAM file.
            sub_bam:  Path to subsampled BAM file.
            bed_dir:  Directory to fetch a random BED file from.
            bed:      Specific BED file to plot.
            out_plt:  Path for the output plot.
        """
        logger.info("Plotter - Initialize")

        self.intervals = Intervals(bed_dir=bed_dir, bed_file=bed)
        self.boundaries = set(
            list(self.intervals.bed["start"]) + list(self.intervals.bed["end"])
        )

        self.in_bam = Loader(bam_path=in_bam)
        self.map_bam = Loader(bam_path=map_bam)
        self.sub_bam = Loader(bam_path=sub_bam)
        self.bams = {"in": self.in_bam, "map": self.map_bam, "out": self.sub_bam}

        self.colormap = {"in": "#648FFF", "map": "#DC267F", "out": "#FFB000"}

        super().check_file_exists(path=out_plt)
        self.out_plt = out_plt

        logger.info("Plotter - Complete initialization")

    def get_pileups(self) -> list:
        """Pileup BAMs for the defined region"""
        logger.info("Plotter - Pileup BAMs")

        pileups = [
            bam.bam.pileup(
                contig=bam.normalize_contig(contig=self.intervals.contig),
                start=self.intervals.start,
                end=self.intervals.end,
            )
            for bam in self.bams.values()
        ]

        logger.info("Plotter - Complete pileup BAMs")
        return pileups

    def get_counts(self, start: int, end: int) -> list:
        """Get read counts from BAMs for the defined region"""
        logger.info("Plotter - Pileup BAMs")

        counts = [
            bam.bam.count(
                contig=bam.normalize_contig(self.intervals.contig),
                start=start,
                end=end,
            )
            for bam in self.bams.values()
        ]

        logger.info("Plotter - Complete pileup BAMs")
        return list(counts)

    def plot(self) -> None:
        """Plot provided BAM file pileups and read counts"""
        logger.info("Plotter - Begin plotting")
        fig, ax_line, ax_bar = self.setup_plot()

        self.add_lineplot(ax=ax_line)
        self.add_boundaries(ax=ax_line)

        self.add_fractions(ax=ax_line)
        self.add_barplot(ax=ax_bar)

        self.add_annotations(ax_line=ax_line, ax_bar=ax_bar)

        logger.info("Plotter - Save plot")
        plt.savefig(self.out_plt, dpi=600)

    @staticmethod
    def setup_plot() -> tuple:
        """Set up the figure and axes for plotting"""
        logger.info("Plotter - Setup plot")
        fig, ax_line = plt.subplots(layout="constrained")
        ax_bar = ax_line.twinx()
        return fig, ax_line, ax_bar

    def add_lineplot(self, ax) -> None:
        """Add line plot representing coverage depth to supplied axis"""
        logger.info("Plotter - Plot pileups")

        for p, l, c in zip(
            self.get_pileups(), self.bams.keys(), self.colormap.values()
        ):
            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )
            ax.plot(
                pileup["coord"], pileup["depth"], label=l, alpha=0.8, color=c, zorder=1
            )

    def add_boundaries(self, ax) -> None:
        """Add vertical lines signifying interval boundaries to supplied axis"""
        logger.info("Plotter - Draw boundaries")

        for b in self.boundaries:
            ax.axvline(
                x=b, linestyle="--", linewidth=1.5, color="gray", alpha=0.5, zorder=1
            )

    def add_annotations(self, ax_line, ax_bar) -> None:
        """Add various text labels to supplied axes"""
        logger.info("Plotter - Add axes, title, legend")

        ax_line.ticklabel_format(useOffset=False, style="plain")
        ax_line.set_title(
            f"Coverage across {self.intervals.contig}:{self.intervals.start}-{self.intervals.end}"
        )
        ax_line.set_xlabel("Chromosomal coordinate")
        ax_line.set_ylabel("Depth of coverage")
        ax_bar.set_ylabel("Read count")
        ax_line.legend(loc="upper right")

    def add_fractions(self, ax):
        """Add fractions to represent distribution values to supplied axis"""
        logger.info("Plotter - Add fraction annotations")

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
        """Add barplot signifying number of reads in each interval"""
        logger.info("Plotter - Add barplot")

        for row in self.intervals.bed.iterrows():
            counts = self.get_counts(start=row[1]["start"], end=row[1]["end"])

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
