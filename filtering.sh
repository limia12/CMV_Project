#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <quality_threshold> <num_cores>"
    exit 1
fi

# Assign arguments to variables
input_dir="$1"
quality_threshold="$2"
num_cores="$3"

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
  echo "Error: The directory $input_dir does not exist."
  exit 1
fi

# Create output directories or warn if they already exist
filtered_dir="$input_dir/filtered_reads"
fastqc_dir="$input_dir/fastqc_results"

if [ -d "$filtered_dir" ]; then
  echo "Warning: $filtered_dir already exists. Existing files will be replaced."
  rm -rf "$filtered_dir"
fi
mkdir -p "$filtered_dir"

if [ -d "$fastqc_dir" ]; then
  echo "Warning: $fastqc_dir already exists. Existing files will be replaced."
  rm -rf "$fastqc_dir"
fi
mkdir -p "$fastqc_dir"

# Paired-end trimming and renaming
echo "Starting paired-end trimming using Trim Galore..."

export quality_threshold filtered_dir

find "$input_dir" -name "*_R1.fastq" | parallel -j $num_cores '
  r1={};
  r2=$(echo {} | sed "s/_R1.fastq/_R2.fastq/");
  base=$(basename {} _R1.fastq);

  trim_galore --quality $quality_threshold --paired --cores 1 \
    --output_dir "$filtered_dir" "$r1" "$r2";

  # Rename output to keep _R1.fastq / _R2.fastq format
  mv "$filtered_dir/${base}_R1_val_1.fq" "$filtered_dir/${base}_R1.fastq";
  mv "$filtered_dir/${base}_R2_val_2.fq" "$filtered_dir/${base}_R2.fastq";
'

echo "Trimming complete. Starting FastQC analysis."

# Run FastQC on all renamed trimmed files
find "$filtered_dir" -name "*.fastq" | parallel -j $num_cores fastqc -t 1 --outdir "$fastqc_dir" {}

echo "FastQC analysis complete. Results saved in $fastqc_dir."
echo "Filtered reads saved in $filtered_dir."
