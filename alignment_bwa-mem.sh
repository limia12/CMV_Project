#!/bin/bash

# Usage:
# alignment_bwa-mem.sh <fastq_dir> <genome_index_dir> <output_root> [--threads <num_cores>]

set -euo pipefail

threads=4

# Parse arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <fastq_dir> <genome_index_dir> <output_root> [--threads <num_cores>]" >&2
    exit 1
fi

fastq_dir=$(realpath "$1")
genome_dir=$(realpath "$2")
output_root=$(realpath "$3")

if [ "${4:-}" == "--threads" ]; then
    if [[ "${5:-}" =~ ^[0-9]+$ ]]; then
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

# Tools check
for tool in bwa samtools; do
    command -v "$tool" >/dev/null || { echo "$tool not found"; exit 1; }
done

bam_dir="$output_root/BAM"
bam_stats_dir="$bam_dir/BAM_Stats"
mkdir -p "$bam_stats_dir"

# Find BWA index prefix
index_prefix=$(find "$genome_dir" -name "*.amb" | head -n 1 | sed 's/\.amb$//')
[ -z "$index_prefix" ] && echo "No .amb found in '$genome_dir'" && exit 1

cd "$fastq_dir" || exit 1

for fq1 in *_R1.fastq; do
    fq2="${fq1/_R1.fastq/_R2.fastq}"
    [ ! -f "$fq2" ] && echo "⚠️  Paired file for $fq1 not found. Skipping." && continue

    sample=$(basename "$fq1" _R1.fastq)
    echo "🔄 Aligning: $sample (threads=$threads)"

    temp_bam="$bam_dir/${sample}.bam"
    sorted_bam="$bam_dir/${sample}_sorted.bam"

    bwa mem -M -t "$threads" "$index_prefix" "$fq1" "$fq2" | samtools view -bS - > "$temp_bam"
    samtools sort "$temp_bam" -o "$sorted_bam"
    samtools index "$sorted_bam"

    samtools flagstat "$sorted_bam" > "$bam_stats_dir/${sample}_flagstat.txt"
    samtools stats "$sorted_bam" > "$bam_stats_dir/${sample}_stats.txt"

    rm -f "$temp_bam"  # 🔥 Clean up intermediate
    echo "✅ Finished: $sample"
done

echo "🎉 All alignments complete. Sorted BAMs in '$bam_dir'."
