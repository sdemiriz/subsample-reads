from subsample_reads.BAMloader import BAMloader
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
        start: str,
        end: str,
        interval_length: str | None,
        interval_count: str | None,
        bed_filename: str,
    ) -> None:
        """
        Constructor for class
        """
        info("Mapper - Initialize BAMcharter")

        self.bam_filename = bam_filename
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)
        self.bed_filename = bed_filename

        if interval_length:
            self.interval_length = int(interval_length)
            self.interval_count = None
        if interval_count:
            self.interval_length = None
            self.interval_count = int(interval_count)

        info("Mapper - Initialize BAMloader")
        self.bam = BAMloader(file=self.bam_filename)

        info("Mapper - Get intervals for BED file")
        self.bed = self.construct_intervals()

        self.total_read_count = self.bam.bam.count(
            contig=self.contig, start=self.start, end=self.end
        )

        info("Mapper - Calculate fraction of reads included in each interval")
        self.bed["fraction"] = [
            self.get_fraction(start=row[1], end=row[2])
            for row in self.bed.itertuples(index=False)
        ]
        self.bed["fraction"] = self.bed["fraction"] / sum(self.bed["fraction"])

        info("Mapper - Write interval data to BED file")
        self.write_bed()

    def get_fraction(self, start: int, end: int) -> float:
        """
        Get number of reads in interval out of all reads in file
        """
        info(
            f"Mapper - Get read count in {start}-{end} interval as fraction of total reads"
        )
        return (
            self.bam.bam.count(contig=self.contig, start=start, end=end)
            / self.total_read_count
        )

    def construct_intervals(self) -> pd.DataFrame:
        """
        Construct BED-formatted DataFrame
        """
        info("Mapper - Form intervals for BED file")
        interval_boundaries = self.get_interval_boundaries()

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
            columns=bed_columns,
        )

    def get_interval_boundaries(self) -> list[int]:
        """
        Divide region based on interval size or count
        """
        region_length = self.end - self.start

        if self.interval_length:
            info("Mapper - Using interval size to subdivide region")

            interval_boundaries = [
                i + self.start for i in range(0, region_length, self.interval_length)
            ]
            if interval_boundaries[-1] != region_length:
                interval_boundaries.append(self.end)

        if self.interval_count:
            info("Mapper - Using interval count to subdivide region")

            interval_size = round(region_length / self.interval_count)
            interval_boundaries = [
                self.start + (i * interval_size)
                for i in range(0, self.interval_count + 1)
            ]
            interval_boundaries[-1] = self.end

        return interval_boundaries

    def write_bed(self) -> None:
        """
        Write output to BED file using filename specified
        """
        info("Mapper - Write BED contents to file")
        self.bed.to_csv(self.bed_filename, sep="\t", index=False, header=False)
