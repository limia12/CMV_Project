#!/bin/bash

set -euo pipefail

# Input directory paths are provided by the GUI
FASTQ_DIR="$1"
BLAST_DB_DIR="$2"
OUT_DIR="$3"
THREADS="$4"

FASTA_DIR="${OUT_DIR}/fasta"
BLAST_OUT_DIR="${OUT_DIR}/blast_output"

mkdir -p "$FASTA_DIR" "$BLAST_OUT_DIR"

# Check dependencies
for cmd in blastn seqtk zcat; do
  if ! command -v $cmd &>/dev/null; then
    echo "❌ Error: '$cmd' not found in PATH. Please install it."
    exit 1
  fi
done

# Loop through FASTQ or FASTQ.GZ files
for fq in "$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fastq.gz; do
  [ -e "$fq" ] || continue
  base=$(basename "$fq")
  base="${base%%.*}"  # strip extensions like .fastq or .fastq.gz
  fasta_file="${FASTA_DIR}/${base}.fasta"

  echo "🔄 Processing sample: $base"

  # Convert FASTQ to FASTA
  if [[ "$fq" == *.gz ]]; then
    zcat "$fq" | seqtk seq -a - > "$fasta_file"
  else
    seqtk seq -a "$fq" > "$fasta_file"
  fi

  # Run BLASTn
  echo "🧬 Running BLAST for $base..."
  blastn -query "$fasta_file" \
         -db "${BLAST_DB_DIR}/viral_nt_db" \
         -evalue 1e-5 \
         -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
         -num_threads "$THREADS" \
         > "${BLAST_OUT_DIR}/${base}_blast.tsv"
done

echo "✅ All BLAST runs complete. Output located in: $BLAST_OUT_DIR"
echo "➡️  Next step: Run the Python summarizer on the BLAST output directory."
