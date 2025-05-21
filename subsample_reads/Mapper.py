from subsample_reads.Loader import Loader
from logging import info
from pathlib import Path
import pandas as pd
import pysam


class Mapper:
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_paths: list[str],
        contig: str,
        start: str,
        end: str,
        interval_length: str | None,
        interval_count: None | str,
        bed_dir: str | None,
    ) -> None:
        """
        Constructor for class
        """
        info("Mapper - Initialize Mapper")

        self.bam_paths = bam_paths
        self.bed_dir = bed_dir
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)

        self.make_bed_dir()
        self.make_bed_filenames()
        self.handle_intervals(
            interval_length=interval_length, interval_count=interval_count
        )
        self.load_bams()
        self.construct_beds()
        self.populate_read_counts()
        self.populate_fractions()

        info("Mapper - Write interval data to BED files")
        self.write_beds()

    def make_bed_dir(self) -> None:
        """
        Make the target directory to place BED files into
        """
        self.bed_dir = Path(self.bed_dir)
        self.bed_dir.mkdir(parents=True, exist_ok=True)

    def make_bed_filenames(self) -> None:
        """
        Generate output BED filenames from  BAM filenames
        """
        self.bed_filenames = [
            self.bed_dir / Path(bam_path).with_suffix(".bed").name
            for bam_path in self.bam_paths
        ]

    def handle_intervals(self, interval_length: int, interval_count: int) -> None:
        """
        Settle whether to use interval length or count when mapping
        """
        if interval_length:
            self.interval_length = int(interval_length)
            self.interval_count = None
        elif interval_count:
            self.interval_length = None
            self.interval_count = int(interval_count)

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

        bed_columns = ["contig", "start", "end", "read_count", "fraction"]

        self.beds = []
        for bam in self.bams:
            intervals = []
            for start, end in zip(interval_boundaries[:-1], interval_boundaries[1:]):
                intervals.append(
                    {
                        bed_columns[0]: self.contig,
                        bed_columns[1]: start,
                        bed_columns[2]: end,
                        bed_columns[3]: -1,
                        bed_columns[4]: -1,
                    }
                )

            bed = pd.DataFrame.from_records(
                intervals,
                columns=bed_columns,
            )

            self.beds.append(bed)

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

    def populate_fractions(self):
        """
        Fill read counts in all BED DataFrames
        """
        info(f"Mapper - Populate read counts in BED files")
        for bed, bam in zip(self.beds, self.bams):
            bed["fraction"] = [
                bam.bam.count(contig=self.contig, start=row[1], end=row[2])
                / sum(bed["read_count"])
                for row in bed.itertuples(index=False)
            ]

    def write_beds(self) -> None:
        """
        Write output to BED file using filename specified
        """
        info("Mapper - Write BED contents to file")
        for bed, path in zip(self.beds, self.bed_filenames):
            bed.to_csv(path, sep="\t", index=False, header=False)
