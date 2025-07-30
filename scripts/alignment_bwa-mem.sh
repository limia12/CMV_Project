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

shopt -s nullglob

matched=false
for fq1 in *_R1.fastq.gz *_R1.fastq *_R1.fq.gz *_R1.fq; do
    matched=true
    fq2="${fq1/_R1./_R2.}"
    [ ! -f "$fq2" ] && echo "⚠️  Paired file for $fq1 not found. Skipping." && continue

    sample=$(basename "$fq1" | sed -E 's/_R1\.(fastq|fq)(\.gz)?$//')
    echo "[$(date +'%F %T')] 🔄 Aligning: $sample (threads=$threads)"

    sorted_bam="$bam_dir/${sample}_sorted.bam"

    # Stream bwa → samtools sort; use samtools threads for speed
    bwa mem -M -t "$threads" "$index_prefix" "$fq1" "$fq2" \
      | samtools sort -@ "$threads" -o "$sorted_bam" -

    samtools index -@ "$threads" "$sorted_bam"
    samtools flagstat "$sorted_bam" > "$bam_stats_dir/${sample}_flagstat.txt"
    samtools stats    "$sorted_bam" > "$bam_stats_dir/${sample}_stats.txt"

    echo "[$(date +'%F %T')] ✅ Finished: $sample"
done

if [ "$matched" = false ]; then
    echo "No FASTQs matched in '$fastq_dir' (expected *_R1.fastq[.gz] or *_R1.fq[.gz])."
fi

echo "🎉 All alignments complete. Sorted BAMs in '$bam_dir'."
