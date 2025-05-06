#!/bin/bash

# Prompt for input directories
read -p "Enter the directory of the FASTQ files: " fastq_dir
read -p "Enter the directory of the genome index files: " genome_dir
read -p "Number of threads to use [default: 4]: " threads
threads=${threads:-4}

# Ask for output directory location
read -p "Enter the directory where you want the BAM output folder to be created: " output_root

# Validate input directories
if [ ! -d "$fastq_dir" ]; then
    echo "Error: FASTQ directory $fastq_dir does not exist."
    exit 1
fi

if [ ! -d "$genome_dir" ]; then
    echo "Error: Genome index directory $genome_dir does not exist."
    exit 1
fi

if [ ! -d "$output_root" ]; then
    echo "Error: Output root directory $output_root does not exist."
    exit 1
fi

# Create/replace BAM and BAM_Stats directories
bam_dir="$output_root/BAM"
bam_stats_dir="$bam_dir/BAM_Stats"

rm -rf "$bam_dir"
mkdir -p "$bam_stats_dir"

# Identify BWA index prefix from .amb file
INDEX=$(find "$genome_dir" -name "*.amb" | sed 's/\.amb$//' | head -n 1)
if [ -z "$INDEX" ]; then
    echo "Error: Could not find a BWA index (.amb) file in $genome_dir."
    exit 1
fi

echo "Using genome index: $INDEX"
cd "$fastq_dir" || exit

# Process all paired FASTQ files
for r1 in *_R1.fastq *_1.fastq; do
    [ -e "$r1" ] || continue

    if [[ "$r1" == *_R1.fastq ]]; then
        base="${r1%_R1.fastq}"
        r2="${base}_R2.fastq"
    elif [[ "$r1" == *_1.fastq ]]; then
        base="${r1%_1.fastq}"
        r2="${base}_2.fastq"
    fi

    if [ -f "$r2" ]; then
        echo "Aligning sample: $base"
        bwa mem -M -k 16 -t "$threads" "$INDEX" "$r1" "$r2" | \
        samtools view -b - | \
        samtools addreplacerg -r "@RG\tID:$base\tSM:$base\tPL:illumina" - > "$bam_dir/${base}.bam"

        # Sort and index BAM
        samtools sort -@ "$threads" -o "$bam_dir/${base}_sorted.bam" "$bam_dir/${base}.bam"
        mv "$bam_dir/${base}_sorted.bam" "$bam_dir/${base}.bam"
        samtools index "$bam_dir/${base}.bam"

        # Generate alignment statistics
        samtools flagstat "$bam_dir/${base}.bam" > "$bam_stats_dir/${base}_flagstat.txt"
        samtools stats "$bam_dir/${base}.bam" > "$bam_stats_dir/${base}_stats.txt"
        samtools idxstats "$bam_dir/${base}.bam" > "$bam_stats_dir/${base}_idxstats.txt"

        echo "Alignment complete for $base"
        echo "Stats stored in: $bam_stats_dir"
    else
        echo "Warning: Paired file not found for $r1 (expected $r2). Skipping."
    fi
done
