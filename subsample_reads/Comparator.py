from logging import info

import pysam
import pandas as pd

from subsample_reads.Loader import Loader


class Comparator:

    def __init__(self, bam_left_path: str, bam_right_path: str, out: str) -> None:
        """ """
        info(
            f"Comparator - Initialize Comparator with {bam_left_path} and {bam_right_path}"
        )

        self.fields_of_interest = [
            "name",
            "ref_name",
            "ref_pos",
            "length",
            "next_ref_name",
            "next_ref_pos",
            # "is_unmapped",
            # "mapping_quality",
            # "is_read1",
            # "is_read2",
            # "is_proper_pair",
            # "mate_is_unmapped",
            # "is_reverse",
            # "mate_is_reverse",
            # "is_secondary",
            # "is_supplementary",
            # "is_duplicate",
            # "flag",
            # "template_length",
            "query_length",
            "query_sequence",
            "query_qualities",
            # "cigarstring",
        ]

        bam_left = Loader(bam_left_path)
        bam_right = Loader(bam_right_path)

        # Get all query_names from the smaller BAM
        bam_right_query_names = set()
        for read in bam_right.fetch():
            bam_right_query_names.add(read.query_name)

        bam_right_query_names_list = list(bam_right_query_names)

        info(
            f"Comparator - {len(bam_right_query_names)} query names in {bam_right_path}"
        )

        # Use hash-based approach for faster matching
        bam_left_df = self.get_matching_reads_hash(bam_left, bam_right_query_names)

        info(
            f"Comparator - Found {len(bam_left_df)}/{len(bam_right_query_names)} reads in {bam_left_path}"
        )

        # Use hash-based approach for right BAM as well for consistency
        bam_right_df = self.get_matching_reads_hash(bam_right, bam_right_query_names)

        info(
            f"Comparator - Found {len(bam_right_df)}/{len(bam_right_query_names)} reads in {bam_right_path}"
        )

        info(f"Comparator - Inner join BAMs")

        merged_df = pd.merge(
            left=bam_left_df,
            right=bam_right_df,
            how="right",
            on="name",
            suffixes=("_l", "_r"),  # type: ignore
        )

        # merged_df["pos_diff"] = int(merged_df["ref_pos_l"]) - int(
        #     merged_df["ref_pos_r"]
        # )
        # merged_df["end_diff"] = merged_df["ref_end_l"] - merged_df["ref_end_r"]

        merged_df.to_csv(out, sep="\t", index=False)

    def get_matching_reads_hash(self, bam: Loader, target_names: set) -> pd.DataFrame:
        """
        Get matching reads using hash-based lookup for O(1) performance
        """
        data = []
        for read in bam.fetch():
            if read.query_name in target_names:
                data.append(self.get_fields_of_interest(read))
        return pd.DataFrame.from_records(data)

    def get_matching_reads(self, bam: Loader) -> pd.DataFrame:
        """
        Get the matching reads from a BAM file
        """
        data = []
        for read in bam.fetch():
            data.append(self.get_fields_of_interest(read))
        return pd.DataFrame.from_records(data)

    def get_fields_of_interest(self, read: pysam.AlignedSegment) -> dict:
        """
        Get the fields of interest from the read
        """
        return {k: v for k, v in read.to_dict().items() if k in self.fields_of_interest}
