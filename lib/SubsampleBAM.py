import pysam
import matplotlib.pyplot as plt
import math


class SubsampleBAM:

    def __init__(
        self,
        in_file: str,
        out_file: str,
        coords: list[tuple[str, int, int, float]],
        seed: int = 42,
        limits: tuple[str, int, int, float] = ("6", 0, 171_115_067, 1.0),
    ):
        """
        Initialize class using input and output BAM files, genomic coordinates, subsetting coefficient, seed
        """
        np.random.normal(loc=1, scale=1, size=12)

        # File IO
        self.in_file = in_file
        self.out_file = out_file

        # Output filename check
        self.out_sorted_file = f"{''.join(self.out_file.split('.')[:-1])}-sorted.bam"
        assert (
            self.out_sorted_file != self.in_file
        ), f"Change {out_file=} to not overwrite {in_file=}"

        # Outer limits for region plotted
        self.limits = limits
        assert (
            len(limits) == 4
        ), f"Number of limits should be exactly 4, not {len(self.limits)}"

        # Define what regions to plot
        self.coords = coords
        for coord in self.coords:
            assert (
                coord[0] == self.limits[0]
            ), f"Supplied chromosome numbers do not match {self.limits[0]=} != {self.coord[0]}"
            assert (
                coord[2] > coord[1]
            ), f"Start {coord[2]} value should be less than end {coord[1]} value"
            assert (
                coord[3] >= 0.0 and coord[3] <= 1.0
            ), f"{coord[3]=} outside the range 0-1 (inclusive)"

        # Define which region to plot using what coefficient
        self.segments = []
        points = (
            [(limits[0], limits[1], limits[1], limits[3])]
            + coords
            + [(limits[0], limits[2], limits[2], limits[3])]
        )
        for i, point in enumerate(points):
            try:
                self.segments.append(
                    (points[i][0], points[i][1] + 1, points[i + 1][1], points[i][3])
                )
            except IndexError as ie:
                self.segments.pop()
                self.segments.append(
                    (points[-1][0], points[i - 1][1] + 1, points[-1][2], points[-1][3])
                )

        assert (
            len(self.segments) == len(self.coords) + 1
        ), f"Number of segments ({len(self.segments)}) != number of coords + 1 ({len(self.coords)})"

        # Set seed
        self.seed = seed

    def run(self, filter_mates: bool = True, plot: bool = False):
        """
        Run the seeded subsampling procedure
        """
        np.random.seed(self.seed)

        with pysam.AlignmentFile(self.in_file, "rb") as in_bam, pysam.AlignmentFile(
            self.out_file, "wb", template=in_bam
        ) as out_bam:

            for seg in self.segments:
                reads = [read for read in in_bam.fetch(seg[0], seg[1], seg[2])]
                print(f"Number of reads in region: {len(reads)}")

                subset_reads = np.random.choice(
                    reads, size=math.ceil(seg[3] * len(reads)), replace=False
                )

                if filter_mates:
                    subset_mates = []
                    nomate_count = 0
                    for read in subset_reads:
                        try:
                            subset_mates.append(in_bam.mate(reads))
                        except:
                            nomate_count += 1

                    print(f"No-mates read count: {nomate_count}\n")

                    for read in subset_reads:
                        if read not in subset_mates:
                            out_bam.write(read)

                else:
                    for read in subset_reads:
                        out_bam.write(read)

        pysam.sort(self.out_file, "-o", self.out_sorted_file)
        pysam.index(self.out_sorted_file)

        if plot:
            self.plot()

    def plot(self):
        """
        Plot in and out BAM file depths
        """
        with pysam.AlignmentFile(self.in_file, "rb") as in_bam, pysam.AlignmentFile(
            self.out_sorted_file, "rb"
        ) as out_bam:

            in_pileup = [
                (p.pos, p.n)
                for p in in_bam.pileup(self.limits[0], self.limits[1], self.limits[2])
            ]
            out_pileup = [
                (p.pos, p.n)
                for p in out_bam.pileup(self.limits[0], self.limits[1], self.limits[2])
            ]

            fig, ax = plt.subplots(layout="constrained")

            ax.plot(
                [p[0] for p in in_pileup],
                [p[1] for p in in_pileup],
                label=f"{self.in_file.split('/')[-1]}",
            )
            ax.plot(
                [p[0] for p in out_pileup],
                [p[1] for p in out_pileup],
                label=f"{self.out_sorted_file.split('/')[-1]}",
            )

            ax.grid(visible=True, linestyle="--", linewidth=1)
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.set_title(
                f"Coverage across chr{self.limits[0]}:{self.limits[1]}-{self.limits[2]}"
            )
            ax.set_xlabel("Chromosomal coordinate")
            ax.set_ylabel("Depth of coverage")
            ax.legend()

            plt.show()


# in_file = "temp-downsample/AH8VC6ADXX-sorted.bam"
# in_file = "results/samples/HG002-hs37d5-300x.bam"
# out_file = "temp-downsample/out.bam"
# coords = [("6", 29_909_037, 29_913_661, 0.1), ("6", 29_913_662, 29_918_286, 0.1)]
# limits = ("6", coords[0][1] - 1000, coords[1][2] + 1000, 1.0)

# test = SubsampleBAM(in_file=in_file, out_file=out_file, coords=coords, limits=limits)
# test.run(filter_mates=True, plot=True)
