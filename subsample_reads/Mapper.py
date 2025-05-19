from subsample_reads.Loader import Loader
from logging import info
from pathlib import Path
import pandas as pd


class Mapper:
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_paths: list[Path],
        contig: str,
        start: int,
        end: int,
        interval_length: int | None,
        interval_count: None | int,
        bed_paths: list[Path],
    ) -> None:
        """
        Constructor for class
        """
        info("Mapper - Initialize Mapper")

        # Receive BAM and BED paths
        self.bam_paths = bam_paths
        self.bed_paths = bed_paths

        # Receive mapping coordinates
        self.contig = contig
        self.start = start
        self.end = end
        self.interval_length = interval_length
        self.interval_count = interval_count

        # Process BAM files into BEDs
        self.load_bams()
        self.construct_beds()
        self.populate_read_counts()

        # Write BEDs to file
        self.write_beds()

        info("Mapper - Complete Mapper")

    def load_bams(self) -> None:
        """
        Initialize Loaders for all supplied BAMs
        """
        info(f"Mapper - Initialize Loaders for supplied BAM files")
        self.bams = [Loader(file=bam_path) for bam_path in self.bam_paths]

    def construct_beds(self) -> pd.DataFrame:
        """
        Construct BED-formatted DataFrame
        """
        info("Mapper - Form intervals for BED file")
        interval_boundaries = self.get_interval_boundaries()

        # For all BAMs, create dataframes with start and end coordinates of each interval
        self.beds = []
        for bam in self.bams:
            intervals = []
            for start, end in zip(interval_boundaries[:-1], interval_boundaries[1:]):
                intervals.append(
                    {
                        "contig": self.contig,
                        "start": start,
                        "end": end,
                        "read_count": -1,
                    }
                )

            bed = pd.DataFrame.from_records(intervals)
            self.beds.append(bed)

    def get_interval_boundaries(self) -> list[int]:
        """
        Divide region based on interval size or count
        """
        region_length = self.end - self.start

        if self.interval_length:
            info("Mapper - Using interval length to subdivide region")

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

    def populate_read_counts(self):
        """
        Fill read counts in all BED DataFrames
        """
        info(f"Mapper - Populate read counts in BED files")
        for bed, bam in zip(self.beds, self.bams):
            bed["read_count"] = [
                bam.bam.count(contig=self.contig, start=row[1], end=row[2])
                for row in bed.itertuples(index=False)
            ]

    def write_beds(self) -> None:
        """
        Write output to BED file using filename specified
        """
        info("Mapper - Write BED contents to file")
        for bed, path in zip(self.beds, self.bed_paths):
            bed.to_csv(path, sep="\t", index=False, header=False)
