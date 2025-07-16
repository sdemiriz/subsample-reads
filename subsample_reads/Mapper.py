from typing import Optional
from pathlib import Path
from math import log
import pandas as pd
import logging
import pysam

from subsample_reads.FileHandler import FileHandler
from subsample_reads.Loader import Loader

logger = logging.getLogger(__name__)


class Mapper(FileHandler):
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_paths: list[str],
        contig: str,
        start: str,
        end: str,
        interval_length: Optional[str] = None,
        interval_count: Optional[str] = None,
        bed_dir: str = "bed/",
    ) -> None:
        """
        Constructor for class
        """
        logger.info("Mapper - Initialize Mapper")

        self.make_bed_dir(path=bed_dir)

        self.bam_paths = bam_paths
        self.bed_paths = self.get_bed_paths(bed_dir=self.bed_dir)

        self.interval_length, self.interval_count = self.setup_intervals(
            interval_length=interval_length,
            interval_count=interval_count,
        )

        self.bams = self.load_bams(bam_paths=self.bam_paths)

        self.contig = contig
        self.beds = self.make_beds(
            contig=self.contig,
            start=start,
            end=end,
            interval_length=self.interval_length,
            interval_count=self.interval_count,
        )

        self.populate_read_counts(contig=self.contig, beds=self.beds, bams=self.bams)
        self.populate_fractions(contig=self.contig, beds=self.beds, bams=self.bams)

        self.write_beds(beds=self.beds, bed_paths=self.bed_paths)

    def make_bed_dir(self, path: str) -> None:
        """
        Make the target directory to place BED files into
        """
        self.bed_dir = Path(path)
        self.bed_dir.mkdir(parents=True, exist_ok=True)

    def get_bed_paths(self, bed_dir: Path) -> list[Path]:
        """
        Generate output BED filenames from  BAM filenames
        """
        return [
            bed_dir / Path(bam_path).with_suffix(".bed").name
            for bam_path in self.bam_paths
        ]

    def setup_intervals(
        self, interval_length: Optional[str], interval_count: Optional[str]
    ) -> tuple[Optional[int], Optional[int]]:
        """
        Settle whether to use interval length or count when mapping
        """
        if interval_length:
            length = int(interval_length)
            count = None
        elif interval_count:
            length = None
            count = int(interval_count)
        else:
            e = "Mapper - No interval length or count provided"
            logger.error(e)
            raise ValueError(e)

        return length, count

    def load_bams(self, bam_paths: list[str]) -> list[Loader]:
        """
        Initialize Loaders for all supplied BAMs
        """
        logger.info(f"Mapper - Initialize Loaders for supplied BAM files")
        return [Loader(bam_path=path) for path in bam_paths]

    def make_beds(
        self,
        contig: str,
        start: str,
        end: str,
        interval_length: Optional[int] = None,
        interval_count: Optional[int] = None,
    ) -> list[pd.DataFrame]:
        """
        Construct BED-formatted DataFrame
        """
        logger.info("Mapper - Form intervals for BED file")
        interval_boundaries = self.get_interval_boundaries(
            start=int(start),
            end=int(end),
            interval_length=interval_length,
            interval_count=interval_count,
        )

        bed_columns = ["contig", "start", "end", "read_count", "fraction"]

        beds = []
        for bam in self.bams:
            intervals = []
            for s, e in zip(interval_boundaries[:-1], interval_boundaries[1:]):
                intervals.append(
                    {
                        bed_columns[0]: contig,
                        bed_columns[1]: s,
                        bed_columns[2]: e,
                        bed_columns[3]: -1,
                        bed_columns[4]: -1,
                    }
                )

            bed = pd.DataFrame.from_records(
                intervals,
                columns=bed_columns,
            )

            beds.append(bed)

        return beds

    def get_interval_boundaries(
        self,
        start: int,
        end: int,
        interval_length: Optional[int],
        interval_count: Optional[int],
    ) -> list[int]:
        """
        Divide region based on interval size or count
        """
        region_length = end - start

        if interval_length:
            logger.info("Mapper - Use interval size to set up intervals")

            interval_boundaries = [
                i + start for i in range(0, region_length, interval_length)
            ]
            if interval_boundaries[-1] != region_length:
                interval_boundaries.append(end)

        elif interval_count:
            logger.info("Mapper - Use interval count to set up intervals")

            interval_size = round(region_length / interval_count)
            interval_boundaries = [
                start + (i * interval_size) for i in range(0, interval_count + 1)
            ]
            interval_boundaries[-1] = end

        return interval_boundaries

    def populate_read_counts(
        self, contig: str, beds: list[pd.DataFrame], bams: list[Loader]
    ) -> None:
        """
        Fill read counts in all BED DataFrames
        """
        logger.info(f"Mapper - Populate read counts in BED files")
        for bed, bam in zip(beds, bams):
            bed["read_count"] = [
                bam.bam.count(contig=contig, start=row[1], end=row[2])
                for row in bed.itertuples(index=False)
            ]

    def populate_fractions(
        self, contig: str, beds: list[pd.DataFrame], bams: list[Loader]
    ) -> None:
        """
        Fill read fractions in all BED DataFrames
        """
        logger.info(f"Mapper - Populate read fractions in BED files")
        for bed, bam in zip(beds, bams):
            bed["fraction"] = [
                bam.bam.count(contig=contig, start=row[1], end=row[2])
                / sum(bed["read_count"])
                for row in bed.itertuples(index=False)
            ]

    def write_beds(self, beds: list[pd.DataFrame], bed_paths: list[Path]) -> None:
        """
        Write output to BED file using filename specified
        """
        logger.info("Mapper - Write BED contents to file")
        for bed, path in zip(beds, bed_paths):
            bed.to_csv(path, sep="\t", index=False, header=False)
            super().check_file_exists(path=str(path))
