from subsample_reads.Intervals import Intervals
from subsample_reads.Plotter import Plotter
from subsample_reads.Loader import Loader
from subsample_reads.Mapper import Mapper

from pathlib import Path
import numpy as np


class Subsampler:

    def __init__(
        self,
        sample_bam_paths: list[str],
        map_bam_paths: list[str],
        bed_dir: str,
        bed_files: list[str] | None,
        bed_count: None | int,
        contig: str,
        start: int,
        end: int,
        interval_length: int | None,
        interval_count: None | int,
        seed: int,
        out_dir: str,
    ):
        # Set main seed for later sampling
        self.main_seed = int(seed)

        # Pathify supplied directories
        self.bed_dir = Path(bed_dir)
        self.out_dir = Path(out_dir)

        # Generate BAM, BED and plot paths
        self.in_bam_paths = [Path(p) for p in sample_bam_paths]
        self.out_bam_paths = [
            self.in_bam_to_out_bam(bam_path=p) for p in sample_bam_paths
        ]

        self.map_bam_paths = [Path(p) for p in map_bam_paths]
        self.bed_paths = self.handle_beds(bed_files=bed_files, bed_count=bed_count)
        self.plt_paths = [self.bam_to_plt(bam_path=p) for p in self.in_bam_paths]

        # Get interval components
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)
        self.handle_intervals(
            interval_length=interval_length, interval_count=interval_count
        )

        # Create BED files for each supplied BAM file
        self.mapper = Mapper(
            bam_paths=self.map_bam_paths,
            contig=self.contig,
            start=self.start,
            end=self.end,
            interval_length=self.interval_length,
            interval_count=self.interval_count,
            bed_paths=self.bed_paths,
        )

        self.intervals = Intervals(bed_paths=self.bed_paths, main_seed=self.main_seed)

        # Run sampling for each input BAM file
        for in_bam_path, out_bam_path, plt_path in zip(
            self.in_bam_paths, self.out_bam_paths, self.plt_paths
        ):
            loader = Loader(path=in_bam_path)
            loader.sample(
                intervals=self.intervals, main_seed=self.main_seed, out_bam=out_bam_path
            )
            loader.close()

            Plotter(
                bam_paths=[in_bam_path] + [out_bam_path],
                intervals=self.intervals,
                out_path=plt_path,
            )

    def in_bam_to_out_bam(self, bam_path: Path) -> Path:
        return self.out_dir / Path(self.add_coords(path=bam_path))

    def bam_to_bed(self, bam_path: Path) -> Path:
        return self.bed_dir / Path(self.add_coords(path=bam_path)).with_suffix(".bed")

    def bam_to_plt(self, bam_path: Path) -> Path:
        return self.out_dir / Path(self.add_coords(path=bam_path)).with_suffix(".png")

    def add_coords(self, path: Path):
        return f"{path.name}-{self.contig}-{self.start}-{self.end}" + path.suffix

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

    def handle_beds(self, bed_files: list[str], bed_count: int):
        if bed_files:
            beds = [self.bed_dir / bed_file for bed_file in bed_files]
        elif bed_count:
            np.random.seed(self.main_seed)
            beds = [
                self.bed_dir
                / np.random.choice(
                    a=list(self.bed_dir.glob("*.bed")), count=bed_count, replace=False
                )
            ]

        return beds
