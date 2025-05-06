#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bam_dir> <out_base>"
    exit 1
fi

# Assign arguments to variables
bam_dir="$1"
out_base="$2"

# Full path for the unaligned reads output folder
unaligned_dir="$out_base/unaligned_reads"

# Create or replace the output directory
if [ -d "$unaligned_dir" ]; then
    echo "Directory '$unaligned_dir' already exists. Replacing it..."
    rm -rf "$unaligned_dir"
fi
mkdir -p "$unaligned_dir"

# Loop over BAM files in the input directory
for bam_file in "$bam_dir"/*.bam; do
    if [ -f "$bam_file" ]; then
        sample_name=$(basename "$bam_file" .bam)
        echo "Processing $sample_name..."

        # Extract unaligned reads (both mates unmapped)
        samtools view -b -f 12 -F 256 "$bam_file" > "$unaligned_dir/${sample_name}_unaligned.bam"

        # Convert unaligned BAM to paired FASTQ files (renamed to R1/R2.fastq)
        samtools fastq \
            -1 "$unaligned_dir/${sample_name}_R1.fastq" \
            -2 "$unaligned_dir/${sample_name}_R2.fastq" \
            -0 /dev/null -s /dev/null -n \
            "$unaligned_dir/${sample_name}_unaligned.bam"

        # Clean up intermediate BAM
        rm "$unaligned_dir/${sample_name}_unaligned.bam"
    else
        echo "No BAM files found in $bam_dir."
    fi
done

echo "All unaligned reads have been extracted to: $unaligned_dir"
