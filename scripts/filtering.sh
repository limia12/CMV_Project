#!/bin/bash

# Usage: filtering.sh <input_dir> <output_dir> <quality_threshold> <num_cores>

# Check for correct number of arguments
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <input_dir> <output_dir> <quality_threshold> <num_cores>"
  exit 1
fi

input_dir="$1"
output_dir="$2"
quality_threshold="$3"
num_cores="$4"

# Validate input directory
if [ ! -d "$input_dir" ]; then
  echo "Error: Input directory '$input_dir' does not exist."
  exit 1
fi

# Prepare output directories
filtered_dir="$output_dir/filtered_reads"
fastqc_dir="$output_dir/fastqc_results"

echo "Setting up directories..."

mkdir -p "$filtered_dir"
mkdir -p "$fastqc_dir"

# Clear old contents if directories exist
rm -rf "$filtered_dir"/*
rm -rf "$fastqc_dir"/*

echo "Starting paired-end trimming using Trim Galore..."
echo "Quality threshold: $quality_threshold | Cores: $num_cores"

export quality_threshold filtered_dir

# Trim reads in parallel
find "$input_dir" -name "*_R1.fastq" | parallel -j "$num_cores" '
  r1={};
  r2=$(echo "$r1" | sed "s/_R1.fastq/_R2.fastq/");
  base=$(basename "$r1" _R1.fastq);

  echo "Processing $base..."

  trim_galore --quality "$quality_threshold" --paired --cores 1 --output_dir "$filtered_dir" "$r1" "$r2"

  mv "$filtered_dir/${base}_R1_val_1.fq" "$filtered_dir/${base}_R1.fastq"
  mv "$filtered_dir/${base}_R2_val_2.fq" "$filtered_dir/${base}_R2.fastq"
'

echo "Trimming complete. Starting FastQC..."

# Run FastQC
find "$filtered_dir" -name "*.fastq" | parallel -j "$num_cores" fastqc -t 1 --outdir "$fastqc_dir" {}

echo "FastQC complete."
echo "Filtered reads saved in: $filtered_dir"
echo "FastQC results saved in: $fastqc_dir"
