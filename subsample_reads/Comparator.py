from subsample_reads.Loader import Loader
from logging import info
import pandas as pd


class Comparator:

    def __init__(self, bam1_path: str, bam2_path: str, out: str) -> None:
        """ """
        info(f"Comparator - Initialize Comparator")

        bam1 = Loader(bam1_path)
        bam2 = Loader(bam2_path)

        # Get all query_names from the smaller BAM
        bam2_query_names = set()
        for read in self.bam3.fetch():
            bam2_query_names.add(read.query_name)
        # If bam2 is smaller, swap the logic

        # Build DataFrame for bam2, filtering by query_name
        bam1_data = []
        for read in self.bam1.fetch():
            if read.query_name in bam2_query_names:
                bam1_data.append(self._extract_read_info(read, self.bam1))
        bam1_df = pd.DataFrame.from_records(bam1_data)

        # Build DataFrame for bam1 (all, or filter as well)
        bam2_data = []
        for read in self.bam2.fetch():
            if read.query_name in bam1_df['query_name'].values:
                bam2_data.append(self._extract_read_info(read, self.bam2))
        bam2_df = pd.DataFrame.from_records(bam2_data)

        info(f"Comparator - Load {bam1_path} and {bam2_path}")

        bam1_df = self.reads_to_df(bam=bam1)
        bam2_df = self.reads_to_df(bam=bam2)

        info(f"Comparator - Inner join BAMs")

        merged_df = pd.merge(
            left=bam1_df,
            right=bam2_df,
            how="inner",
            on="query_name",
            suffixes=["_1", "_2"],
        )

        merged_df["start_diff"] = merged_df["ref_start_1"] - merged_df["ref_start_2"]
        merged_df["end_diff"] = merged_df["ref_end_1"] - merged_df["ref_end_2"]

        info(f"Comparator - Write to {out}")
        merged_df.to_csv(out, sep="\t", index=False)

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
