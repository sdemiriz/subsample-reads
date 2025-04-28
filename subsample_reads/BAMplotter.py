from subsample_reads.BAMloader import BAMloader
from subsample_reads.Intervals import Intervals
import matplotlib.pyplot as plt
from logging import info
import pandas as pd
import pysam


class BAMplotter:

    def __init__(
        self,
        bam_files: list[str],
        bed_file: str,
        out: str,
    ):
        """
        Constructor for plotting utility
        """
        info(f"Initialize BAMplotter with {bam_files=}, {bed_file=}, {out=}")

        self.bed_file = bed_file
        self.bed = Intervals(self.bed_file)

        self.bam_files = bam_files
        self.get_depths()

        self.out = out
        self.plot()

        info(f"Complete BAMplotter")

    def get_depths(self):
        """
        Pileup BAMs for the defined region
        """
        info("Pileup BAMs")

        self.depth_files = [bam.split(".")[0] + ".csv" for bam in self.bam_files]
        for bam, depth in zip(self.bam_files, self.depth_files):
            pysam.depth(bam, "-r", "6", "-o", depth)

        info("Complete pileup BAMs")

    def plot(self):
        """
        Plot provided BAM file pileups
        """
        info(f"Begin plotting")
        fig, ax = plt.subplots(layout="constrained")

        contig, start, end = self.bed.get_limits()

        for depth in self.depth_files:

            depths = pd.read_csv(depth, sep="\t", header=None, usecols=[1, 2])
            depths = depths[::1000]

            ax.plot(
                depths[1],
                depths[2],
                label=f"{depth.split('.')[-2]}",
            )

        ax.set_yscale("log")
        ax.grid(visible=True, linestyle="--", linewidth=1)
        # ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across {contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        info(f"Complete plotting")

        info(f"Save plot")
        plt.savefig(self.out)
