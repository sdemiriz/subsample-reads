import pysam

# Make a simple example BAM file with 100 reads, 100bp each, on chr1
# The reads are spaced 100bp apart, and the BAM file is 1000bp long
# The BAM file is indexed after creation

# Header creation
header = pysam.AlignmentHeader.from_dict(
    {
        "HD": {"VN": "1.5"},
        "PG": [],
        "SQ": [{"SN": "chr1", "LN": 248956422}],
        "RG": [],
    }
)
out_file = "example.bam"
out_bam = pysam.AlignmentFile(out_file, "wb", header=header)

# Read creation
read_count = 0
for interval in range(10):
    interval_start = interval * 100 + 1
    interval_end = interval_start + 100

    for i in range(10):
        read = {
            "name": f"READ_{read_count:04d}",
            "ref_name": "chr1",
            "flag": "0",
            "ref_pos": str(interval_start),
            "next_ref_pos": str(interval_start + 1),
            "next_ref_name": "chr1",
            "cigar": "100M",
            "seq": "A" * 100,
            "qual": "?" * 100,
            "tags": [],
            "map_quality": "60",
            "length": "1000",
        }

        new_read = pysam.AlignedSegment.from_dict(read, header=header)
        out_bam.write(new_read)

        read_count += 1

out_bam.close()
pysam.index(out_file)
