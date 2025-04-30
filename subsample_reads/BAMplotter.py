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
        out: str | None,
    ) -> None:
        """
        Constructor for plotting utility
        """
        info(f"Initialize BAMplotter with {bam_files=}, {bed_file=}, {out=}")

        self.bed_file = bed_file
        self.bed = Intervals(file=self.bed_file)

        self.bam_files = bam_files

        self.out = out
        self.plot()

        info(f"Complete BAMplotter")

    def get_pileups(self) -> list:
        """
        Pileup BAMs for the defined region
        """
        info("Pileup BAMs")

        contig, start, end = self.bed.get_limits()
        bams = [BAMloader(file=bam) for bam in self.bam_files]
        pileups = [
            bam.bam.pileup(contig=bam.handle_contig_name(contig), start=start, end=end)
            for bam in bams
        ]

        info("Complete pileup BAMs")
        return pileups

    def plot(self) -> None:
        """
        Plot provided BAM file pileups
        """
        info(f"Begin plotting")
        fig, ax = plt.subplots(layout="constrained")

        contig, start, end = self.bed.get_limits()
        pileups = self.get_pileups()

        for p in pileups:

            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )

            ax.plot(
                pileup["coord"],
                pileup["depth"],
                label=f"{self.bed_file.split('.')[-2]}",
            )

        ax.grid(visible=True, linestyle="--", linewidth=1)
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across {contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        info(f"Complete plotting")

        info(f"Save plot")
        plt.savefig(self.out)
