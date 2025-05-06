#!/bin/bash

# Usage:
# alignment_bwa-mem.sh <fastq_dir> <genome_index_dir> <output_root> [--threads <num_cores>]

set -euo pipefail

# Default number of threads
threads=4

# Parse arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <fastq_dir> <genome_index_dir> <output_root> [--threads <num_cores>]" >&2
    exit 1
fi

fastq_dir=$(realpath "$1")
genome_dir=$(realpath "$2")
output_root=$(realpath "$3")

# Optional: parse --threads
if [ "${4:-}" == "--threads" ]; then
    if [ -n "${5:-}" ] && [[ "$5" =~ ^[0-9]+$ ]]; then
        threads="$5"
    else
        echo "Error: --threads requires a numeric argument." >&2
        exit 1
    fi
fi

# Validate input directories
for dir in "$fastq_dir" "$genome_dir" "$output_root"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' does not exist." >&2
        exit 1
    fi
done

# Check required tools
for tool in bwa samtools; do
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: '$tool' not found in PATH." >&2
        exit 1
    fi
done

# Setup output dirs
bam_dir="$output_root/BAM"
bam_stats_dir="$bam_dir/BAM_Stats"
rm -rf "$bam_dir"
mkdir -p "$bam_stats_dir"

# Get BWA index prefix
index_prefix=$(find "$genome_dir" -name "*.amb" | head -n 1 | sed 's/\.amb$//')
if [ -z "$index_prefix" ]; then
    echo "Error: No .amb index found in '$genome_dir'." >&2
    exit 1
fi

echo "🧬 Using BWA index: $index_prefix"
cd "$fastq_dir" || exit 1

# Align paired FASTQs
for fq1 in *_R1.fastq; do
    fq2="${fq1/_R1.fastq/_R2.fastq}"
    if [ ! -f "$fq2" ]; then
        echo "⚠️  Skipping: Paired file for '$fq1' not found."
        continue
    fi

    sample=$(basename "$fq1" _R1.fastq)
    echo "🔄 Aligning: $sample (threads=$threads)"

    bwa mem -M -t "$threads" "$index_prefix" "$fq1" "$fq2" | \
        samtools view -bS - > "$bam_dir/${sample}.bam"

    samtools sort "$bam_dir/${sample}.bam" -o "$bam_dir/${sample}_sorted.bam"
    samtools index "$bam_dir/${sample}_sorted.bam"

    samtools flagstat "$bam_dir/${sample}_sorted.bam" > "$bam_stats_dir/${sample}_flagstat.txt"
    samtools stats "$bam_dir/${sample}_sorted.bam" > "$bam_stats_dir/${sample}_stats.txt"

    echo "✅ Finished: $sample"
done

echo "🎉 All alignments complete. BAMs in '$bam_dir', stats in '$bam_stats_dir'."
