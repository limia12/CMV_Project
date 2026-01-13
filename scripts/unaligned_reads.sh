#!/bin/bash
# unaligned_reads.sh
#
# Purpose:
#   For each BAM in <bam_dir>, extract read pairs where BOTH mates are unmapped
#   and write them as paired FASTQ files (R1 and R2).
#
# Usage:
#   unaligned_reads.sh <bam_dir> <out_base> [threads]
#
# Inputs:
#   bam_dir   = directory containing input *.bam files
#   out_base  = base output directory (writes to <out_base>/unaligned_reads)
#   threads   = optional samtools thread count (default: 1)
#
# Outputs:
#   <out_base>/unaligned_reads/<sample>_R1.fastq
#   <out_base>/unaligned_reads/<sample>_R2.fastq

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Inputs ------------------------------ #
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <bam_dir> <out_base> [threads]"
    exit 1
fi

bam_dir="$1"          # Directory containing input BAM files.
out_base="$2"         # Base output directory.
threads="${3:-1}"     # Threads for samtools (default: 1).

# -------------------------- Output directory ------------------------ #
unaligned_dir="$out_base/unaligned_reads"

# Remove output directory if it already exists to avoid mixing old and new files.
if [ -d "$unaligned_dir" ]; then
    echo "Directory '$unaligned_dir' already exists. Replacing it..."
    rm -rf "$unaligned_dir"
fi

# Create output directory (and parents if needed).
mkdir -p "$unaligned_dir"

# If no BAMs match, do not treat the pattern as a literal string.
shopt -s nullglob

# ------------------------- Per-sample loop -------------------------- #
for bam_file in "$bam_dir"/*.bam; do
    [ -f "$bam_file" ] || continue

    # Sample name is the BAM filename without the path and ".bam" extension.
    sample_name=$(basename "$bam_file" .bam)
    echo "[INFO] Processing sample: $sample_name"

    # 1) Extract read pairs where both reads are unmapped.
    #
    # Flag logic used by samtools view:
    #   -f 12  requires BOTH flags:
    #          4 = read is unmapped
    #          8 = mate is unmapped
    #   -F 256 excludes secondary alignments.
    #
    # Output is a temporary BAM that only contains unmapped pairs.
    samtools view -@ "$threads" -b -f 12 -F 256 "$bam_file" \
        > "$unaligned_dir/${sample_name}_unaligned.bam"

    # 2) Convert the temporary BAM into paired FASTQ files.
    #
    # -1 writes read 1 to R1 FASTQ.
    # -2 writes read 2 to R2 FASTQ.
    # -0 /dev/null discards reads that do not fit into R1/R2 outputs.
    # -s /dev/null discards singletons (mates missing).
    # -n keeps read names consistent for pairing.
    samtools fastq -@ "$threads" \
        -1 "$unaligned_dir/${sample_name}_R1.fastq" \
        -2 "$unaligned_dir/${sample_name}_R2.fastq" \
        -0 /dev/null -s /dev/null -n \
        "$unaligned_dir/${sample_name}_unaligned.bam"

    # 3) Delete the temporary BAM to save disk space.
    rm "$unaligned_dir/${sample_name}_unaligned.bam"
done

# ------------------------------ Done -------------------------------- #
echo "[INFO] All unaligned reads have been extracted to: $unaligned_dir"
