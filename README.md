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
python count_deletions.py input.bam 250 > deletion_counts.tsv
```

## find_deletions.py
Returns a table of all deletions found in a BAM file with their size and location.

Usage: 
```
pip install pysam
python find_deletions.py input.bam > deletions.tsv
```
