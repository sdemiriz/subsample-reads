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
            "query_name",
            "ref_name",
            "ref_start",
            "ref_end",
            "next_ref_name",
            "next_ref_start",
            "is_unmapped",
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
            # "query_length",
            # "query_sequence",
            # "query_qualities",
            # "cigarstring",
        ]

        bam_left = Loader(bam_left_path)
        bam_right = Loader(bam_right_path)

        # Get all query_names from the smaller BAM - use more efficient set operations
        bam_right_query_names = set()
        for read in bam_right.fetch():
            bam_right_query_names.add(read.query_name)
        # Convert to list only if needed for other operations
        bam_right_query_names_list = list(bam_right_query_names)

        info(
            f"Comparator - {len(bam_right_query_names)} query names in {bam_right_path}"
        )

        bam_left_df = self.get_matching_reads(bam_left)

        info(
            f"Comparator - Found {len(bam_left_df)}/{len(bam_right_query_names)} reads in {bam_left_path}"
        )

        bam_right_df = self.get_matching_reads(bam_right)

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

        info(f"Comparator - Write to {out}")
        merged_df.to_csv(out, sep="\t", index=False)

    def get_matching_reads(self, bam: Loader) -> pd.DataFrame:
        """
        Get the matching reads from the two BAM files
        """
        data = []
        for read in bam.fetch():
            data.append(self.get_fields_of_interest(read))
        return pd.DataFrame.from_records(data)

    def get_fields_of_interest(self, read: pysam.AlignedSegment) -> dict[str, Any]:
        """
        Get the fields of interest from the read
        """
        return {k: v for k, v in read.to_dict().items() if k in self.fields_of_interest}

    def reads_to_df(self, bam: Loader) -> pd.DataFrame:
        """
        Read contents from loaded BAM file into Pandas DataFrame
        """
        read_data = []

        # Get all read contents
        for read in bam.fetch():
            read_info = {
                "query_name": read.query_name,
                "ref_name": (
                    bam.get_reference_name(read.reference_id)
                    if read.reference_id != -1
                    else None
                ),
                "ref_start": read.reference_start,
                "ref_end": read.reference_end,
                "next_ref_name": (
                    bam.get_reference_name(read.next_reference_id)
                    if read.next_reference_id != -1
                    else None
                ),
                "next_ref_start": read.next_reference_start,
                "is_unmapped": read.is_unmapped,
                # "mapping_quality": read.mapping_quality,
                # "is_read1": read.is_read1,
                # "is_read2": read.is_read2,
                # "is_proper_pair": read.is_proper_pair,
                # "mate_is_unmapped": read.mate_is_unmapped,
                # "is_reverse": read.is_reverse,
                # "mate_is_reverse": read.mate_is_reverse,
                # "is_secondary": read.is_secondary,
                # "is_supplementary": read.is_supplementary,
                # "is_duplicate": read.is_duplicate,
                # "flag": read.flag,
                # "template_length": read.template_length,
                # "query_length": read.query_length,
                # "query_sequence": read.query_sequence,
                # "query_qualities": read.query_qualities,
                # "cigarstring": read.cigarstring,
            }

        # Get tag values
        for tag, value in read.get_tags():
            read_info[tag] = value

        # Append to list before emitting DataFrame
        read_data.append(read_info)

        return pd.DataFrame.from_dict(read_data)
