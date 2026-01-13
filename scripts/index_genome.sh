#!/bin/bash
# index_genome.sh
#
# Purpose:
#   Prepare a reference genome for common NGS tools by creating:
#     - BWA index files (for read mapping)
#     - FASTA index (.fai) using samtools (for random access)
#     - Sequence dictionary (.dict) using Picard (often required by GATK)
#
# Usage:
#   index_genome.sh <genome_fasta_or_fasta.gz>
#
# Input:
#   A reference genome FASTA file, either plain (.fa/.fasta) or gzipped (.gz).
#
# Output:
#   Index files are saved alongside the original FASTA file, using the same
#   base name (even if the input was gzipped).

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Inputs ------------------------------ #
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <genome_fasta_or_fasta.gz>"
    exit 1
fi

genome_fasta="$1"   # Path to the reference FASTA (plain or gzipped).

# Check the FASTA exists.
if [ ! -f "$genome_fasta" ]; then
    echo "Error: The genome FASTA file $genome_fasta does not exist."
    exit 1
fi

# ------------------------- Output naming ---------------------------- #
# Work out where to save files (same directory as the input FASTA).
dest_dir="$(dirname "$(realpath "$genome_fasta")")"

# If the input ends with .gz, strip it so outputs match the FASTA base name.
base_no_gz="${genome_fasta%.gz}"

# Prefix used for output files (.dict, .fai, and BWA index files).
dest_prefix="${dest_dir}/$(basename "$base_no_gz")"

# --------------------- Handle gzipped FASTA ------------------------- #
# BWA requires an uncompressed FASTA. If the input is gzipped, decompress
# to a temporary file and use that for indexing.
work_fa="$base_no_gz"   # The FASTA we will actually index.
tmpdir=""               # Only used if we need to decompress.

if [[ "$genome_fasta" == *.gz ]]; then
    tmpdir="$(mktemp -d)"
    work_fa="$tmpdir/$(basename "$base_no_gz")"
    echo "[INFO] Decompressing $genome_fasta -> $work_fa"
    gzip -cd "$genome_fasta" > "$work_fa"
fi

# ------------------------- Build indexes ---------------------------- #
# 1) BWA index: creates files needed for read alignment.
echo "[INFO] Running BWA index..."
bwa index "$work_fa"

# 2) samtools faidx: creates a .fai file used for fast access to FASTA.
echo "[INFO] Creating FASTA index (.fai) with samtools..."
samtools faidx "$work_fa"

# 3) Picard sequence dictionary: creates a .dict file often required by GATK.
echo "[INFO] Creating sequence dictionary (.dict) with Picard..."
picard CreateSequenceDictionary R="$work_fa" O="${dest_prefix}.dict"

# --------------------- Move outputs into place ---------------------- #
# If we used a temporary FASTA, the index files are created in tmpdir.
# Move them next to the original FASTA and use the original base name.
for ext in amb ann bwt pac sa; do
    src="${work_fa}.${ext}"
    if [ -f "$src" ]; then
        mv -f "$src" "${dest_prefix}.${ext}"
    fi
done

# Move the FASTA index file.
if [ -f "${work_fa}.fai" ]; then
    mv -f "${work_fa}.fai" "${dest_prefix}.fai"
fi

# ------------------------------ Cleanup ----------------------------- #
# Remove temporary folder if we created one.
if [ -n "$tmpdir" ]; then
    rm -rf "$tmpdir"
fi

echo "[INFO] Genome indexing complete. Files saved next to the FASTA."
