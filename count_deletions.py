#!/usr/bin/env python3
import pysam
import sys
from collections import Counter

def count_deletions(bam_file, min_read_length=0):
    # Dictionary to store counts of deletion sizes
    deletion_sizes = Counter()

    # Open BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
                
            # Skip reads shorter than minimum length
            if read.query_length < min_read_length:
                continue

            # Iterate over CIGAR operations
            for op, length in read.cigartuples or []:
                # CIGAR code 2 means deletion from reference
                if op == 2:
                    deletion_sizes[length] += 1

    return deletion_sizes


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.bam> <min_read_length>")
        sys.exit(1)

    bam_file = sys.argv[1]
    min_read_length = int(sys.argv[2])
    deletions = count_deletions(bam_file, min_read_length)

    bam_name = bam_file.split('/')[-1]
    print(f"DeletionSize\t{bam_name}")
    for size in sorted(deletions):
        print(f"{size}\t{deletions[size]}")
