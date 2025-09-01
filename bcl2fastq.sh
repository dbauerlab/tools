#!/bin/bash

# Script for running bcl2fastq conversion

# Slurm job parameters
#SBATCH --job-name=bcl2fastq_conversion
#SBATCH --output=bcl2fastq_conversion.out
#SBATCH --error=bcl2fastq_conversion.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load necessary modules
ml bcl2fastq2/2.20.0-GCC-12.2.0

# Define input and output directories
INPUT_DIR="/path/to/run/folder"
OUTPUT_DIR="/path/to/fastq/output"
SAMPLE_SHEET="/path/to/SampleSheet.csv"

# Run bcl2fastq conversion
bcl2fastq -r 1 \
-p 4 \
-w 1 \
--runfolder-dir ${INPUT_DIR} \
--output-dir ${OUTPUT_DIR} \
--sample-sheet ${SAMPLE_SHEET} \
--no-lane-splitting
