#!/bin/bash

# Define reference path
REF_DIR="/home/stp/Bioinformatics/Limia/CMV_Project/References/Simulated_Reads_Ref"
OUT_DIR="/home/stp/Bioinformatics/Limia/CMV_Project/Simulated_Reads"

# Create output directory
mkdir -p "$OUT_DIR"

# Set common ART parameters
ART="/usr/local/bin/art_illumina"  # Update this if needed
READ_LEN=150
FRAG_LEN=300
STD_DEV=10
SEQUENCER="HS25"  # Illumina HiSeq 2500

# Simulate reads and rename to _R1.fastq / _R2.fastq
simulate_reads() {
    local fasta=$1
    local coverage=$2
    local prefix=$3

    # Ensure prefix ends with underscore
    [[ "$prefix" != *_ ]] && prefix="${prefix}_"

    $ART -ss $SEQUENCER \
         -i "$fasta" \
         -p \
         -l $READ_LEN \
         -f $coverage \
         -m $FRAG_LEN \
         -s $STD_DEV \
         -o "${OUT_DIR}/${prefix}"

    # Rename output files
    mv "${OUT_DIR}/${prefix}1.fq" "${OUT_DIR}/${prefix}R1.fastq"
    mv "${OUT_DIR}/${prefix}2.fq" "${OUT_DIR}/${prefix}R2.fastq"
}

echo "Simulating CMV only @100x..."
simulate_reads "$REF_DIR/hcmv_genome.fasta" 100 "CMV_100x"

echo "Simulating chr21 only @100x..."
simulate_reads "$REF_DIR/chr21.fasta" 100 "chr21_100x"

# Simulate chr21 100x + CMV at various coverages
for cov in 1 5 10 25 50; do
    echo "Simulating chr21 @100x + CMV @${cov}x..."
    simulate_reads "$REF_DIR/chr21.fasta" 100 "chr21_100x_CMV_${cov}x_chr21"
    simulate_reads "$REF_DIR/hcmv_genome.fasta" $cov "chr21_100x_CMV_${cov}x_CMV"
done

# Add CMV @100x + chr21 @100x (equal mix)
echo "Simulating chr21 @100x + CMV @100x..."
simulate_reads "$REF_DIR/chr21.fasta" 100 "chr21_100x_CMV_100x_chr21"
simulate_reads "$REF_DIR/hcmv_genome.fasta" 100 "chr21_100x_CMV_100x_CMV"

# Simulate chr21 @100x + each herpesvirus @100x
for virus in HHV-6A HHV-6B HHV-7 HHV-8 EBV; do
    fasta="$REF_DIR/${virus}.fasta"
    echo "Simulating chr21 @100x + ${virus} @100x..."
    simulate_reads "$REF_DIR/chr21.fasta" 100 "chr21_${virus}_chr21"
    simulate_reads "$fasta" 100 "chr21_${virus}_${virus}"
done

echo "✅ All simulations completed. FASTQ files saved in: $OUT_DIR"
