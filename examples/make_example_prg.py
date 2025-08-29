import pysam

# Make a simple example BAM file with 100 reads each on chr6 and PRG_96
# The reads are spaced 100bp apart, and the BAM file region is 1000bp long
# The BAM file is indexed after creation

# Header creation
header = pysam.AlignmentHeader.from_dict(
    {
        "HD": {"VN": "1.5"},
        "PG": [],
        "SQ": [{"SN": "chr6", "LN": 170805979}, {"SN": "PRG_96", "LN": 3517}],
        "RG": [],
    }
)
out_file = "example-prg.bam"
out_bam = pysam.AlignmentFile(out_file, "wb", header=header)

# Read creation on chr6, HLA-A gene
hla_a = ("chr6", 29941260, 29949572)
read_count = 0
for interval in range(10):
    interval_start = hla_a[1] + interval * 100 + 1
    interval_end = interval_start + 100

    for i in range(10):
        read = {
            "name": f"CHR6_READ_{read_count:03d}",
            "ref_name": "chr6",
            "flag": "0",
            "ref_pos": str(interval_start),
            "next_ref_pos": str(interval_start + 1),
            "next_ref_name": "chr6",
            "cigar": "100M",
            "seq": "A" * 100,
            "qual": "?" * 100,
            "tags": [],
            "map_quality": "60",
            "length": "1000",
        }

        read = pysam.AlignedSegment.from_dict(read, header=header)
        out_bam.write(read)

        read_count += 1

# Create reads on PRG_96, which maps back to HLA-A
read_count = 0
for interval in range(10):
    interval_start = interval * 100 + 1
    interval_end = interval_start + 100

    for i in range(10):
        read = {
            "name": f"PRG_READ_{read_count:03d}",
            "ref_name": "PRG_96",
            "flag": "0",
            "ref_pos": str(interval_start),
            "next_ref_pos": str(interval_start + 1),
            "next_ref_name": "PRG_96",
            "cigar": "100M",
            "seq": "A" * 100,
            "qual": "?" * 100,
            "tags": [],
            "map_quality": "60",
            "length": "1000",
        }

        read = pysam.AlignedSegment.from_dict(read, header=header)
        out_bam.write(read)

        read_count += 1

out_bam.close()
pysam.index(out_file)
