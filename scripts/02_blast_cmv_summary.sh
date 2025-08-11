#!/usr/bin/env bash
set -euo pipefail

# Arguments:
#   $1 = FASTQ_DIR
#   $2 = BLAST_DB_DIR
#   $3 = OUT_DIR  (e.g., /.../output) — single directory; TSV goes to blast_output/, CSVs go here
#   $4 = THREADS

FASTQ_DIR="$1"
BLAST_DB_DIR="$2"
OUT_DIR="$3"
THREADS="$4"

FASTA_DIR="${OUT_DIR}/fasta"
BLAST_OUT_DIR="${OUT_DIR}/blast_output"

mkdir -p "$FASTA_DIR" "$BLAST_OUT_DIR"

# Check dependencies
for cmd in blastn seqtk zcat python3; do
  if ! command -v "$cmd" &>/dev/null; then
    echo "❌ Error: '$cmd' not found. Please install it."
    exit 1
  fi
done

# Run BLAST for each FASTQ/FASTQ.GZ
for fq in "$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fastq.gz; do
  [ -e "$fq" ] || continue
  base=$(basename "$fq")
  base="${base%%.*}"  # strip extensions like .fastq.gz
  fasta_file="${FASTA_DIR}/${base}.fasta"
  echo "🔄 Processing sample: $base"

  if [[ "$fq" == *.gz ]]; then
    zcat "$fq" | seqtk seq -a - > "$fasta_file"
  else
    seqtk seq -a "$fq" > "$fasta_file"
  fi

  echo "🧬 Running BLAST for ${base}..."
  blastn \
    -query "$fasta_file" \
    -db "${BLAST_DB_DIR}/viral_nt_db" \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -num_threads "$THREADS" \
    > "${BLAST_OUT_DIR}/${base}_blast.tsv"
done

echo "✅ All BLAST runs complete. TSVs are in: ${BLAST_OUT_DIR}"
echo "➡️ Summarizing BLAST hits (top5 + heatmap CSVs) into: ${OUT_DIR}"

# Summarize, writing per-sample top5 and heatmap + unified heatmap into OUT_DIR
python3 summarize_blast_hits.py "$BLAST_OUT_DIR" "$OUT_DIR"

echo "✅ BLAST summary (CSV) written to: ${OUT_DIR}"
