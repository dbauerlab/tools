#!/usr/bin/env python3
import pysam
import sys

def analyze_alignments(bam_file, min_read_length=0):
    # Open BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Print header
        print("ReadLength\tChromosome\tAlignStart\tAlignEnd\tMaxDelSize\tDelStart\tDelEnd")
        
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
                
            # Skip reads shorter than minimum length
            if read.query_length < min_read_length:
                continue
            
            # Get basic alignment info
            chrom = bam.get_reference_name(read.reference_id)
            align_start = read.reference_start + 1  # Convert to 1-based
            align_end = read.reference_end
            read_length = read.query_length
            
            # Find the largest deletion in the read
            max_del_size = 0
            max_del_start = 0
            max_del_end = 0
            ref_pos = read.reference_start
            
            # Iterate through CIGAR operations to find deletions
            for op, length in read.cigartuples or []:
                if op == 0 or op == 7 or op == 8:  # Match/Mismatch/Equal
                    ref_pos += length
                elif op == 2:  # Deletion
                    if length > max_del_size:
                        max_del_size = length
                        max_del_start = ref_pos + 1  # Convert to 1-based
                        max_del_end = ref_pos + length
                    ref_pos += length
                elif op == 3:  # N (skip)
                    ref_pos += length
            
            # Output the information
            print(f"{read_length}\t{chrom}\t{align_start}\t{align_end}\t{max_del_size}\t{max_del_start}\t{max_del_end}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.bam> <min_read_length>")
        sys.exit(1)

    bam_file = sys.argv[1]
    min_read_length = int(sys.argv[2])
    
    analyze_alignments(bam_file, min_read_length)
