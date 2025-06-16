#!/bin/bash

# Usage:
# extract_unaligned_reads.sh <bam_dir> <out_base>

# Check that exactly two arguments are provided: input BAM directory and output base path
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bam_dir> <out_base>"
    exit 1
fi

# Assign the input arguments to named variables for readability
bam_dir="$1"           # Directory containing BAM files
out_base="$2"          # Base output path (e.g., project folder)

# Define the full path to the directory where unaligned reads will be saved
unaligned_dir="$out_base/unaligned_reads"

# If the output directory already exists, warn the user and delete it to avoid old results mixing in
if [ -d "$unaligned_dir" ]; then
    echo "Directory '$unaligned_dir' already exists. Replacing it..."
    rm -rf "$unaligned_dir"
fi
# Create a fresh output directory
mkdir -p "$unaligned_dir"

# Loop through all BAM files in the input directory
for bam_file in "$bam_dir"/*.bam; do
    # Ensure it's a valid file (this avoids errors if the directory is empty)
    if [ -f "$bam_file" ]; then
        # Get the filename without the .bam extension to use as sample name
        sample_name=$(basename "$bam_file" .bam)
        echo "Processing $sample_name..."

        # Extract paired reads where BOTH mates are unmapped
        # -f 12 = both reads unmapped, -F 256 = exclude secondary alignments
        samtools view -b -f 12 -F 256 "$bam_file" > "$unaligned_dir/${sample_name}_unaligned.bam"

        # Convert the unaligned BAM file to paired FASTQ files
        # -1: output file for read 1
        # -2: output file for read 2
        # -0: discard reads with only one end unmapped (we want only fully unaligned pairs)
        # -s: discard singleton reads (where mate is missing)
        # -n: ensure proper read name formatting
        samtools fastq \
            -1 "$unaligned_dir/${sample_name}_R1.fastq" \
            -2 "$unaligned_dir/${sample_name}_R2.fastq" \
            -0 /dev/null -s /dev/null -n \
            "$unaligned_dir/${sample_name}_unaligned.bam"

        # Delete the intermediate BAM file to save space
        rm "$unaligned_dir/${sample_name}_unaligned.bam"
    else
        # If no .bam files were found, display a message (this check runs each iteration)
        echo "No BAM files found in $bam_dir."
    fi
done

# Final confirmation message
echo "All unaligned reads have been extracted to: $unaligned_dir"
