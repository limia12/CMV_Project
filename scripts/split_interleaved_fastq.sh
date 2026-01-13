#!/bin/bash
# split_interleaved_fastq.sh
#
# Purpose:
#   Split interleaved FASTQ files into two separate files:
#     - *_R1.fastq (read 1)
#     - *_R2.fastq (read 2)
#
# Notes:
#   - This assumes the input FASTQ is truly interleaved:
#       R1 record (4 lines), then R2 record (4 lines), repeating.
#   - This script processes all "*.fastq" files under <input_dir>.
#
# Usage:
#   split_interleaved_fastq.sh <input_dir>
#
# Output:
#   Creates: <input_dir>/split_interleaved/
#   Writes: <base>_R1.fastq and <base>_R2.fastq for each input FASTQ.

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Inputs ------------------------------ #
# Require exactly one argument: the folder containing interleaved FASTQs.
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_dir>"
    exit 1
fi

input_dir="$1"   # Directory to search for *.fastq files.

# -------------------------- Output directory ------------------------ #
# Create an output folder inside the input directory.
output_dir="$input_dir/split_interleaved"
mkdir -p "$output_dir"

# ---------------------------- Main logic ---------------------------- #
# Split one interleaved FASTQ into *_R1.fastq and *_R2.fastq.
split_fastq() {
    fastq_file="$1"   # Path to the interleaved FASTQ.
    output_dir="$2"   # Where to write the split files.

    # Remove the ".fastq" extension to get a base name for outputs.
    base_name=$(basename "$fastq_file" .fastq)

    # Output filenames.
    r1_file="$output_dir/${base_name}_R1.fastq"
    r2_file="$output_dir/${base_name}_R2.fastq"

    # Create/overwrite output files to start clean.
    : > "$r1_file"
    : > "$r2_file"

    # A FASTQ record is 4 lines:
    #   1) header
    #   2) sequence
    #   3) plus line
    #   4) quality string
    #
    # Interleaved format is:
    #   R1 record (4 lines), then R2 record (4 lines), repeated.
    while read -r l1 && read -r l2 && read -r l3 && read -r l4; do
        # Write the first 4 lines as the next R1 record.
        printf "%s\n%s\n%s\n%s\n" "$l1" "$l2" "$l3" "$l4" >> "$r1_file"

        # Read the next 4 lines as the matching R2 record.
        read -r l1 && read -r l2 && read -r l3 && read -r l4
        printf "%s\n%s\n%s\n%s\n" "$l1" "$l2" "$l3" "$l4" >> "$r2_file"
    done < "$fastq_file"

    echo "[INFO] Split $(basename "$fastq_file") -> $(basename "$r1_file"), $(basename "$r2_file")"
}

# Export the function so it can be used in the subshells started by xargs.
export -f split_fastq

# --------------------------- Batch running -------------------------- #
# Find all "*.fastq" files and run split_fastq on each one.
#
# xargs options:
#   -n 1  = pass one FASTQ path per command
#   -P 16 = run up to 16 files in parallel (adjust if needed)
#   -I {} = replace {} with the FASTQ path in the command
find "$input_dir" -name "*.fastq" \
    | xargs -n 1 -P 16 -I {} bash -c 'split_fastq "$@"' _ {} "$output_dir"
