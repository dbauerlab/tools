# tools
Repo for general and random bioinformatic tools/scripts/functions

## check_fasta_gtf.sh
Checks whether sequence names in FASTA and GTF match 

Usage: 
```
./check_fasta_gtf.sh genome.fasta[.gz] annotations.gtf[.gz]
```

## count_deletions.py
Counts the occurences of deletions of all sizes in a BAM file. Skips reads less than a given length.

Usage: 
```
pip install pysam
python count_deletions.py input.bam 250 > deletion_counts.tsv
```

## find_deletions.py
Returns a table of all deletions found in a BAM file with their size and location. Skips reads less than a given length.

Usage: 
```
pip install pysam
python find_deletions.py input.bam 250 > deletions.tsv
```

## read_alignments.py
For each mapped read that meets the length threshold, it outputs:

- Read length
- Chromosome name
- Alignment start position (1-based)
- Alignment end position
- Size of the largest deletion in the read
- Start position of the largest deletion
- End position of the largest deletion

All positions are 1-based for consistency with standard genomic coordinates. If a read has no deletions, the MaxDelSize will be 0 and DelStart/DelEnd will also be 0. Reads shorter than the specified minimum length are skipped. Unmapped reads are skipped.

Usage:
```
python read_alignments.py input.bam 100 > report.tsv
```
