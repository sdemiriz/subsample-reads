from subsample_reads.BAMloader import BAMloader
from subsample_reads.Intervals import Intervals
import matplotlib.pyplot as plt
from logging import info


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
        self.bams = [BAMloader(bam) for bam in bam_files]

        self.out = out
        self.plot(out=self.out)

        info(f"Complete BAMplotter")

    def pileup_bams(self, contig: str, start: int, end: int):
        """
        Pileup BAMs for the defined region
        """
        info("Pileup BAMs")
        pileups = [
            bam.bam.pileup(contig=contig, start=start, stop=end) for bam in self.bams
        ]

        info("Complete pileup BAMs")
        return pileups

    def plot(self, out):
        """
        Plot provided BAM file pileups
        """
        info(f"Begin plotting")
        fig, ax = plt.subplots(layout="constrained")

        contig, start, end = self.bed.get_limits()
        pileups = self.pileup_bams(contig=contig, start=start, end=end)

        for pileup in pileups:
            ax.plot(
                [p[0] for p in pileup],
                [p[1] for p in pileup],
                label=f"{self.bam.file.split('/')[-1]}",
            )

        ax.grid(visible=True, linestyle="--", linewidth=1)
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title(f"Coverage across {contig}:{start}-{end}")
        ax.set_xlabel("Chromosomal coordinate")
        ax.set_ylabel("Depth of coverage")
        ax.legend()

        info(f"Complete plotting")

        info(f"Save plot")
        plt.savefig(out)
