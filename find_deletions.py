#!/usr/bin/env python3
import pysam
import sys

def find_deletions(bam_file, min_read_length=0):
    # Open BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
                
            # Skip reads shorter than minimum length
            if read.query_length < min_read_length:
                continue

            ref_pos = read.reference_start  # Genomic start position of read

            # Iterate through CIGAR operations
            for op, length in read.cigartuples or []:
                if op == 0 or op == 7 or op == 8:
                    # Match/Mismatch/EQ — advance reference position
                    ref_pos += length
                elif op == 2:
                    # Deletion from reference
                    chrom = bam.get_reference_name(read.reference_id)
                    start = ref_pos
                    end = ref_pos + length
                    yield chrom, start, end, length
                    ref_pos += length
                elif op in (3, 6, 5):
                    # Skips / soft / hard clips: advance ref if needed
                    if op == 3:  # N skip
                        ref_pos += length
                # Other ops (insertions, padding) don't advance reference for deletions

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.bam> <min_read_length>")
        sys.exit(1)

    bam_file = sys.argv[1]
    min_read_length = int(sys.argv[2])

    print("Chromosome\tStart\tEnd\tDeletionSize")
    for chrom, start, end, size in find_deletions(bam_file, min_read_length):
        print(f"{chrom}\t{start+1}\t{end}\t{size}")
