from subsample_reads.BAMloader import BAMloader
import pysam
import pandas as pd
from logging import info


class Mapper:
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_filename: str,
        contig: str,
        start: int,
        end: int,
        interval_length: int | None,
        interval_count: int | None,
        bed_filename: str,
    ) -> None:
        """
        Constructor for class
        """
        info("Initialize BAMcharter")

        self.bam_filename = bam_filename
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)

        self.interval_length = interval_length
        self.interval_count = interval_count
        self.bed_filename = bed_filename

        info("Initialize BAMloader")
        self.bam = BAMloader(file=self.bam_filename)

        info("Get intervals for BED file")
        self.bed = self.construct_intervals()
        self.contig = self.bed["contig"][0]

        info("Calculate fraction of reads included in each interval")
        self.bed["fraction"] = [
            self.get_fraction(begin=row[1], end=row[2])
            for row in self.bed.itertuples(index=False)
        ]

        info(f"Write interval data to BED file")
        self.write_bed()

    def get_num_reads(self):
        """
        Get number of reads in contig
        """
        info("Get number of reads in contig")
        return self.bam.bam.count(contig=self.contig, start=self.start, stop=self.end)

    def get_reads_in_region(self):
        """
        Get read iterator for reads in specified region
        """
        yield from self.bam.bam.fetch(
            contig=self.contig, start=self.start, stop=self.end
        )

    def get_fraction(self, begin: int, end: int):
        """
        Get number of reads in interval out of all reads in file
        """
        info("Get number of reads in interval as fraction of total reads")
        return (
            self.bam.bam.count(contig=self.contig, start=begin, end=end)
            / self.get_num_reads()
        )

    def construct_intervals(self) -> pd.DataFrame:
        """
        Construct BED-formatted DataFrame
        """
        info("Form intervals for BED file")
        interval_boundaries = self.divide_contig()

        bed_columns = ["contig", "start", "end", "fraction"]

        df_precursor = []
        for i in range(len(interval_boundaries)):
            if i == 0:
                pass
            else:
                df_precursor.append(
                    {
                        bed_columns[0]: self.contig,
                        bed_columns[1]: interval_boundaries[i - 1],
                        bed_columns[2]: interval_boundaries[i],
                        bed_columns[3]: -1,
                    }
                )

        return pd.DataFrame.from_records(
            df_precursor,
            columns=["contig", "begin", "end", "fraction"],
        )

    def divide_contig(self):
        """
        Divide full region absed on interval size or interval count
        """
        reference_length = self.bam.bam.get_reference_length(self.contig)

        if self.interval_length:
            info("Using interval size to subdivide region")

            interval_boundaries = [
                i for i in range(0, reference_length, self.interval_length)
            ]
            if interval_boundaries[-1] != reference_length:
                interval_boundaries.append(reference_length)

        if self.interval_count:
            info("Using interval count to subdivide region")

            interval_size = round(reference_length / self.interval_count)
            interval_boundaries = [
                i * interval_size for i in range(0, self.interval_count + 1)
            ]
            interval_boundaries[-1] = reference_length

        return interval_boundaries

    def write_bed(self):
        """
        Write output to BED file using filename specified
        """
        info("Write BED contents to file")
        self.bed.to_csv(self.bed_filename, sep="\t", index=False, header=False)
