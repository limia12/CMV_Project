#!/bin/bash

# Usage:
# trim_and_qc.sh <input_dir> <quality_threshold> <num_cores>

# Check if exactly three arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <quality_threshold> <num_cores>"
    exit 1
fi

# Assign input arguments to meaningful variable names
input_dir="$1"             # Directory with raw FASTQ files
quality_threshold="$2"     # Quality threshold for trimming (e.g., 20)
num_cores="$3"             # Number of CPU cores to use with GNU Parallel

# Check that the input directory exists
if [ ! -d "$input_dir" ]; then
  echo "Error: The directory $input_dir does not exist."
  exit 1
fi

# Define output directories relative to input directory
filtered_dir="$input_dir/filtered_reads"     # Where trimmed FASTQs will be saved
fastqc_dir="$input_dir/fastqc_results"       # Where FastQC results will be saved

# If filtered_reads dir already exists, warn and delete it to avoid mixing results
if [ -d "$filtered_dir" ]; then
  echo "Warning: $filtered_dir already exists. Existing files will be replaced."
  rm -rf "$filtered_dir"
fi
mkdir -p "$filtered_dir"   # Create new filtered_reads directory

# If fastqc_results dir already exists, warn and delete it to avoid conflicts
if [ -d "$fastqc_dir" ]; then
  echo "Warning: $fastqc_dir already exists. Existing files will be replaced."
  rm -rf "$fastqc_dir"
fi
mkdir -p "$fastqc_dir"     # Create new fastqc_results directory

# Notify user that trimming is starting
echo "Starting paired-end trimming using Trim Galore..."

# Export variables so they’re visible to each GNU parallel job
export quality_threshold filtered_dir

# For each *_R1.fastq file found:
# - Determine corresponding *_R2.fastq file
# - Extract sample base name
# - Run trim_galore with paired-end mode on that pair
# - Rename the trimmed output files to standard _R1.fastq / _R2.fastq naming
find "$input_dir" -name "*_R1.fastq" | parallel -j $num_cores '
  r1={};  # Save current R1 filename
  r2=$(echo {} | sed "s/_R1.fastq/_R2.fastq/");  # Infer matching R2 file
  base=$(basename {} _R1.fastq);  # Extract sample ID from filename

  # Run Trim Galore with quality trimming and paired-end mode
  trim_galore --quality $quality_threshold --paired --cores 1 \
    --output_dir "$filtered_dir" "$r1" "$r2";

  # Rename the output back to the expected _R1/_R2.fastq names
  mv "$filtered_dir/${base}_R1_val_1.fq" "$filtered_dir/${base}_R1.fastq";
  mv "$filtered_dir/${base}_R2_val_2.fq" "$filtered_dir/${base}_R2.fastq";
'

# Let user know trimming is complete and FastQC is starting
echo "Trimming complete. Starting FastQC analysis."

# Run FastQC on all trimmed FASTQ files using GNU Parallel
# One thread per job, results go into fastqc_results directory
find "$filtered_dir" -name "*.fastq" | parallel -j $num_cores fastqc -t 1 --outdir "$fastqc_dir" {}

# Final messages to confirm success and where outputs are saved
echo "FastQC analysis complete. Results saved in $fastqc_dir."
echo "Filtered reads saved in $filtered_dir."
