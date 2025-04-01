import pysam


class BAMloader:

    def __init__(self, in_file: str, c: str, start: int, end: int):

        self.in_file = in_file
        self.chr = c
        self.start = start
        self.end = end

        self.load_bam()
        self.reads = self.bam.fetch(contig=self.chr, start=self.start, end=self.end)

    def load_bam(self):
        self.bam = pysam.AlignmentFile(self.in_file, "rb")
