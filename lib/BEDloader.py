import pandas as pd
from intervaltree import Interval, IntervalTree


class BEDloader:

    def __init__(self, in_file: str):

        self.in_file = in_file
        self.load_bed()
        self.create_intervaltree()

        self.tree = IntervalTree()

    def load_bed(self):

        self.bed = pd.read_csv(self.in_file, sep="\t", header=0)
        self.bed.columns = ["chr", "start", "end", "fraction"]

    def create_intervaltree(self):
        for row in self.bed.iterrows():
            self.tree[row["start"], row["end"]] = row["fraction"]
