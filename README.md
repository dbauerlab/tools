# tools
Repo for general and random bioinformatic tools/scripts/functions

## check_fasta_gtf.sh
Checks whether sequence names in FASTA and GTF match 

Usage: 
```
./check_fasta_gtf.sh genome.fasta[.gz] annotations.gtf[.gz]
```

## count_deletions.py
Counts the occurences of deletions of all sizes in a BAM file.

Usage: 
```
pip install pysam
python count_deletions.py input.bam > deletion_counts.tsv
```
