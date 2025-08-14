from logging import info

import pysam
import tempfile
import pandas as pd

from subsample_reads.Loader import Loader


class Comparator:

    def __init__(self, bam_left_path: str, bam_right_path: str, out: str) -> None:
        """ """
        info(f"Comparator - Initialize Comparator")

        self.fields_of_interest = [
            "name",
            "ref_name",
            "ref_pos",
            "length",
            "next_ref_name",
            "next_ref_pos",
        ]

        bam_left = Loader(bam_left_path)
        bam_right = Loader(bam_right_path)

        # Get all query_names from the smaller BAM (right)
        bam_right_query_names = set()
        for read in bam_right.fetch():
            bam_right_query_names.add(read.query_name)

        with tempfile.NamedTemporaryFile(mode="w+") as tmpfile:

            info(f"Comparator - Create temporary file for read names")

            for line in bam_right_query_names:
                tmpfile.write(f"{line}\n")
            tmpfile.flush()

            bam_right_str = pysam.view(str(bam_right_path), "-N", tmpfile.name)
            bam_left_df = self.get_reads(bam=bam_left, read_names_str=bam_right_str)
            bam_right_df = self.get_reads(bam=bam_right, read_names_str=bam_right_str)

        info(f"Comparator - Inner join BAMs")

        suffixes = ("_l", "_r")
        merged_df = pd.merge(
            left=bam_left_df,
            right=bam_right_df,
            how="right",
            on="name",
            suffixes=suffixes,
        )

        column_order = ["name"]
        for col in [c for c in self.fields_of_interest if c != "name"]:
            for suf in suffixes:
                column_order.append(f"{col}{suf}")

        merged_df = merged_df[column_order]
        merged_df.to_csv(out, sep="\t", index=False)

    def get_reads(self, bam: Loader, read_names_str: str) -> pd.DataFrame:

        info(f"Comparator - Get reads from BAM")
        reads = []
        for line in read_names_str.splitlines():
            reads.append(
                pysam.AlignedSegment.fromstring(line, bam.bam.header).to_dict()
            )

        df = pd.DataFrame.from_records(reads)
        df = df[self.fields_of_interest]

        return df
