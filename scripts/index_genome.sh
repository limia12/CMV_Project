#!/bin/bash

# Usage:
# index_genome.sh <genome_fasta>

# Ensure exactly one argument is provided (path to the genome FASTA file)
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <genome_fasta>"
    exit 1
fi

# Assign the argument to a variable for clarity
genome_fasta="$1"

# Check if the specified FASTA file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: The genome FASTA file $genome_fasta does not exist."
  exit 1
fi

# Run BWA to generate genome index files (e.g., .amb, .ann, .bwt, etc.)
# These are needed for alignment with BWA later
echo "Running BWA index on the genome..."
bwa index "$genome_fasta"

# Use Samtools to create a .fai index for the FASTA file
# This helps tools quickly access regions in the reference
echo "Generating FASTA index with Samtools..."
samtools faidx "$genome_fasta"

# Use Picard to create a .dict file (sequence dictionary)
# Required by tools like GATK for reading the genome
echo "Creating sequence dictionary using Picard..."
picard CreateSequenceDictionary R="$genome_fasta" O="${genome_fasta}.dict"

# Final message to confirm indexing is complete
echo "Genome indexing complete. Files saved in the same directory as the FASTA file."
