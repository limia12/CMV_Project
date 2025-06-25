#!/bin/bash

# Usage:
# unaligned_reads.sh <bam_dir> <out_base> [threads]

# Check that at least 2 arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <bam_dir> <out_base> [threads]"
    exit 1
fi

bam_dir="$1"
out_base="$2"
threads="${3:-1}"  # Optional third arg for threads (default: 1)

unaligned_dir="$out_base/unaligned_reads"

if [ -d "$unaligned_dir" ]; then
    echo "Directory '$unaligned_dir' already exists. Replacing it..."
    rm -rf "$unaligned_dir"
fi

mkdir -p "$unaligned_dir"

for bam_file in "$bam_dir"/*.bam; do
    if [ -f "$bam_file" ]; then
        sample_name=$(basename "$bam_file" .bam)
        echo "Processing $sample_name..."

        samtools view -@ "$threads" -b -f 12 -F 256 "$bam_file" > "$unaligned_dir/${sample_name}_unaligned.bam"

        samtools fastq -@ "$threads" \
            -1 "$unaligned_dir/${sample_name}_R1.fastq" \
            -2 "$unaligned_dir/${sample_name}_R2.fastq" \
            -0 /dev/null -s /dev/null -n \
            "$unaligned_dir/${sample_name}_unaligned.bam"

        rm "$unaligned_dir/${sample_name}_unaligned.bam"
    fi
done

echo "All unaligned reads have been extracted to: $unaligned_dir"
