from subsample_reads.BAMloader import BAMloader
import pysam
import pandas as pd
from logging import info


class BAMcharter:
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_file: str,
        bed_file: str,
        contig: str,
        window_size: int,
        window_count: int,
    ):
        """"""
        info(
            f"Initialize BAMcharter with {bam_file=}, {bed_file=}, {contig=}, {window_size=}, {window_count=}"
        )
        self.bam_file = bam_file
        self.bed_file = bed_file

        self.window_size = window_size
        self.window_count = window_count

        info(f"Initialize BAMloader using BAM file {self.bam_file}")
        self.bam = BAMloader(file=self.bam_file)

        info(f"Calculate total number of reads in BAM file")
        self.total_read_count = self.get_num_reads()

        info(f"Generate intervals for BED file using window size")
        self.form_bed_intervals()
        self.contig = self.bed["contig"][0]

        info(f"Calculate fraction of reads included in each interval")
        fractions = []
        for row in self.bed.itertuples(index=False):
            fractions.append(
                self.get_reads_fraction(contig=row[0], begin=row[1], end=row[2])
            )
        self.bed["fraction"] = fractions

        info(f"Write interval data to BED file")
        self.write_bed()

    def get_reads_fraction(self, contig: str, begin: int, end: int):
        """ """
        info(
            f"Calculate number of reads in interval as fraction of total reads in BAM file"
        )
        return (
            self.bam.bam.count(contig=contig, start=begin, end=end)
            / self.total_read_count
        )

    def form_bed_intervals(self):
        """ """
        info(f"Form intervals for BED file")
        interval_boundaries = self.divide_contig()

        to_dataframe = []
        for i, interval in enumerate(interval_boundaries):
            if i == 0:
                pass
            else:
                d = {
                    "contig": self.contig,
                    "begin": interval_boundaries[i - 1],
                    "end": interval_boundaries[i],
                    "fraction": -1,
                }
                to_dataframe.append(d)

        self.bed = pd.DataFrame.from_records(
            to_dataframe,
            columns=["contig", "begin", "end", "fraction"],
        )

    def divide_contig(self):
        """ """
        reference_length = self.bam.bam.get_reference_length(self.contig)

        if self.window_size:
            self.window_size = int(self.window_size)
            info(f"Using {self.window_size=}")

            interval_boundaries = [
                i for i in range(0, reference_length, self.window_size)
            ]
            if interval_boundaries[-1] != reference_length:
                interval_boundaries.append(reference_length)

        if self.window_count:
            self.window_count = int(self.window_count)
            info(f"Using {self.window_count=}")

            interval_size = round(reference_length / self.window_count)
            interval_boundaries = [
                i * interval_size for i in range(0, self.window_count + 1)
            ]
            interval_boundaries[-1] = reference_length

        return interval_boundaries

    def get_num_reads(self):
        """ """
        info(f"Get number of reads in contig {self.contig}")
        return self.bam.bam.count(self.contig)

    def write_bed(self):
        """
        Write output to BED file using filename specified
        """
        info(f"Write BED file contents to {self.bed_file}")
        self.bed.to_csv(self.bed_file, sep="\t", index=False, header=False)
