#!/usr/bin/env python3
import pysam
import sys
from collections import Counter

def count_deletions(bam_file):
    # Dictionary to store counts of deletion sizes
    deletion_sizes = Counter()

    # Open BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            # Iterate over CIGAR operations
            for op, length in read.cigartuples or []:
                # CIGAR code 2 means deletion from reference
                if op == 2:
                    deletion_sizes[length] += 1

    return deletion_sizes


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.bam>")
        sys.exit(1)

    bam_file = sys.argv[1]
    deletions = count_deletions(bam_file)

    bam_name = bam_file.split('/')[-1]
    print(f"DeletionSize\t{bam_name}")
    for size in sorted(deletions):
        print(f"{size}\t{deletions[size]}")
