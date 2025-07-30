#!/bin/bash

# Paths (update if needed)
REF_DIR="/home/stp/Bioinformatics/Limia/CMV_Project/References/Simulated_Reads_Ref"
OUT_DIR="/home/stp/Bioinformatics/Limia/CMV_Project/Simulated_Reads"
ART="/usr/local/bin/art_illumina"

# Common ART settings
READ_LEN=150
FRAG_LEN=300
STD_DEV=10
SEQUENCER="HS25"

# Output prefix
PREFIX="chr21_NonCMV_Herpes_100x"

# Create temp dir for intermediate FASTQs
TMP_DIR="${OUT_DIR}/tmp_${PREFIX}"
mkdir -p "$TMP_DIR"

# Function to simulate and rename
simulate_to_temp() {
    local fasta=$1
    local coverage=$2
    local label=$3

    $ART -ss $SEQUENCER \
         -i "$fasta" \
         -p \
         -l $READ_LEN \
         -f $coverage \
         -m $FRAG_LEN \
         -s $STD_DEV \
         -o "${TMP_DIR}/${label}_"

    # Rename to standard R1/R2 for merge later
    mv "${TMP_DIR}/${label}_1.fq" "${TMP_DIR}/${label}_R1.fastq"
    mv "${TMP_DIR}/${label}_2.fq" "${TMP_DIR}/${label}_R2.fastq"
}

echo "🧬 Simulating chr21 @100x..."
simulate_to_temp "$REF_DIR/chr21.fasta" 100 "chr21"

# Herpesviruses (excluding HCMV)
VIRUSES=("HHV-6A" "HHV-6B" "HHV-7" "HHV-8" "EBV")

for virus in "${VIRUSES[@]}"; do
    echo "🧬 Simulating ${virus} @100x..."
    simulate_to_temp "$REF_DIR/${virus}.fasta" 100 "$virus"
done

echo "🧪 Merging all simulated R1 reads..."
cat "${TMP_DIR}"/*_R1.fastq > "${OUT_DIR}/${PREFIX}_R1.fastq"

echo "🧪 Merging all simulated R2 reads..."
cat "${TMP_DIR}"/*_R2.fastq > "${OUT_DIR}/${PREFIX}_R2.fastq"

echo "🧹 Cleaning up temp files..."
rm -r "$TMP_DIR"

echo "✅ Final paired FASTQ files (No HCMV):"
echo "   ${OUT_DIR}/${PREFIX}_R1.fastq"
echo "   ${OUT_DIR}/${PREFIX}_R2.fastq"
