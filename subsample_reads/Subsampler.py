from subsample_reads.Loader import Loader
from subsample_reads.Mapper import Mapper
from subsample_reads.Plotter import Plotter
from subsample_reads.Intervals import Intervals

import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pysam


class Subsampler:

    def __init__(
        self,
        sample_bam_paths: list[str],
        map_bam_paths: list[str],
        bed_dir: str,
        bed_list: list[str] | None,
        bed_count: None | int,
        contig: str,
        start: int,
        end: int,
        interval_length: int | None,
        interval_count: None | int,
        seed: int,
        out_dir: str,
    ):

        # Pathify supplied directories
        self.bed_dir = Path(bed_dir)
        self.out_dir = Path(out_dir)

        # Generate BAM, BED and plot paths
        self.in_bam_paths = [Path(p) for p in sample_bam_paths]
        self.out_bam_paths = [
            self.in_bam_to_out_bam(bam_path=p) for p in sample_bam_paths
        ]

        self.map_bam_paths = [Path(p) for p in map_bam_paths]
        self.bed_paths = [self.bam_to_bed(bam_path=p) for p in self.in_bam_paths]
        self.plt_paths = [self.bam_to_plt(bam_path=p) for p in self.in_bam_paths]

        # Set main seed for later sampling
        self.main_seed = int(seed)

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

        # Run sampling for each input BAM file
        for in_bam_path, out_bam_path, plt_path in zip(
            self.in_bam_paths, self.out_bam_paths, self.plt_paths
        ):
            loader = Loader(path=in_bam_path)
            loader.sample(
                bed_paths=self.bed_paths, main_seed=self.main_seed, out_bam=out_bam_path
            )
            loader.close()

            Plotter(
                bam_paths=[in_bam_path] + [out_bam_path],
                bed_paths=self.bed_paths,
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
