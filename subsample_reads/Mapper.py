from subsample_reads.Loader import Loader
from logging import info
import pandas as pd
import pysam


class Mapper:
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_filenames: list[str],
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

        self.bam_filenames = bam_filenames
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

        info("Mapper - Initialize Loader")
        self.bams = [Loader(file=bam_file) for bam_file in self.bam_filenames]

        info("Mapper - Get intervals for BED file")
        self.bed = self.construct_intervals()

        for b in self.bams:

            total_read_count = self.get_read_count(
                bam=b, start=self.start, end=self.end
            )
            header_f = f"{b.file}_fraction"
            header_r = f"{b.file}_read_count"

            info("Mapper - Calculate fraction of reads included in each interval")
            self.bed[header_f] = [
                self.get_fraction(
                    bam=b, start=row[1], end=row[2], total_read_count=total_read_count
                )
                for row in self.bed.itertuples(index=False)
            ]
            self.bed[header_f] = self.bed[header_f] / sum(self.bed[header_f])

            self.bed[header_r] = [
                self.get_read_count(bam=b, start=row[1], end=row[2])
                for row in self.bed.itertuples(index=False)
            ]

        self.bed["fraction"] = self.bed[self.filenames_with_suffix("_fraction")].mean(
            axis=1
        )

        self.bed["read_count"] = (
            self.bed[self.filenames_with_suffix("_read_count")].mean(axis=1).astype(int)
        )

        self.bed.drop(
            columns=self.filenames_with_suffix("_fraction")
            + self.filenames_with_suffix("_read_count"),
            inplace=True,
        )

        info("Mapper - Write interval data to BED file")
        self.write_bed()

    def filenames_with_suffix(self, suffix: str):
        return [f + suffix for f in self.bam_filenames]

    def get_fraction(
        self, bam: pysam.AlignmentFile, start: int, end: int, total_read_count: int
    ) -> float:
        """
        Get number of reads in interval out of all reads in file
        """
        info(
            f"Mapper - Get read count in {start}-{end} interval as fraction of total reads"
        )
        return self.get_read_count(bam=bam, start=start, end=end) / total_read_count

    def get_read_count(self, bam: pysam.AlignmentFile, start: int, end: int) -> int:
        """
        Get number of reads in interval out of all reads in file
        """
        info(
            f"Mapper - Get read count in {start}-{end} interval as fraction of total reads"
        )
        return bam.bam.count(contig=self.contig, start=start, end=end)

    def construct_intervals(self) -> pd.DataFrame:
        """
        Construct BED-formatted DataFrame
        """
        info("Mapper - Form intervals for BED file")
        interval_boundaries = self.get_interval_boundaries()

        bed_columns = ["contig", "start", "end", "fraction", "read_count"]

        to_df = []
        for i in range(len(interval_boundaries)):
            if i == 0:
                pass
            else:
                to_df.append(
                    {
                        bed_columns[0]: self.contig,
                        bed_columns[1]: interval_boundaries[i - 1],
                        bed_columns[2]: interval_boundaries[i],
                        bed_columns[3]: -1,
                        bed_columns[4]: -1,
                    }
                )

        return pd.DataFrame.from_records(
            to_df,
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
