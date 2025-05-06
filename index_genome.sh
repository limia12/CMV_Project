#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <genome_fasta>"
    exit 1
fi

# Assign arguments to variables
genome_fasta="$1"

# Check if the genome FASTA file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: The genome FASTA file $genome_fasta does not exist."
  exit 1
fi

# Step 1: Index the genome using BWA
echo "Running BWA index on the genome..."
bwa index "$genome_fasta"

# Step 2: Generate the FASTA index using Samtools
echo "Generating FASTA index with Samtools..."
samtools faidx "$genome_fasta"

# Step 3: Generate the sequence dictionary using Picard
echo "Creating sequence dictionary using Picard..."
picard CreateSequenceDictionary R="$genome_fasta" O="${genome_fasta}.dict"

# Step 4: Confirm completion
echo "Genome indexing complete. Files saved in the same directory as the FASTA file."
