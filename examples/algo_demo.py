import pysam
import random

random.seed(42)

# Make a simple example BAM file with 100 reads, 100bp each, on chr1
# The reads are spaced randomly, and the BAM file is 1000bp long
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
out_file = "examples/algo-demo.bam"
out_bam = pysam.AlignmentFile(out_file, "wb", header=header)
reads = []


# Read creation
read_length = 150
interval_n_start, interval_n_end = 1000, 2000
interval_n_plus_1_start, interval_n_plus_1_end = 2000, 3000


def get_random_start_coords():
    coords_n = [
        random.randint(interval_n_start, interval_n_end - read_length - 1)
        for _ in range(9)
    ]
    coords_n_plus_1 = [
        random.randint(interval_n_plus_1_start, interval_n_plus_1_end - read_length - 1)
        for _ in range(9)
    ]
    for coord in sorted(coords_n + coords_n_plus_1):
        yield coord


start_coords = get_random_start_coords()

read_template = {
    "ref_name": "chr1",
    "flag": "0",
    "next_ref_name": "chr1",
    "cigar": f"{read_length}M",
    "seq": "A" * read_length,
    "qual": "?" * read_length,
    "tags": [],
    "map_quality": "60",
    "length": "1000",
}

# Interval N reads
for i in range(1, 10):

    interval_start = next(start_coords)
    interval_end = interval_start + read_length

    read = read_template.copy()
    read["name"] = f"r{i}"
    read["ref_pos"] = str(interval_start)
    read["next_ref_pos"] = str(interval_start + 1)

    reads.append(read)

# Overhang read
read = read_template.copy()
read["name"] = f"r10"
read["ref_pos"] = str(int(interval_n_end - read_length / 2))
read["next_ref_pos"] = str(int(interval_n_end - read_length / 2 + 1))

reads.append(read)

# Interval N+1 reads
for i in range(11, 20):

    interval_start = next(start_coords)
    interval_end = interval_start + read_length

    read = read_template.copy()
    read["name"] = f"r{i}"
    read["ref_pos"] = str(interval_start)
    read["next_ref_pos"] = str(interval_start + 1)

    reads.append(read)


for read in reads:
    out_bam.write(pysam.AlignedSegment.from_dict(read, header=header))

out_bam.close()
pysam.index(out_file)
