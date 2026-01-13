#!/usr/bin/env bash
# 02_blast_cmv_summary.sh
#
# What this script does:
#   1) Takes FASTQ files (either in a folder, or inside an archive file).
#   2) Converts each FASTQ to FASTA (because BLAST uses FASTA input here).
#   3) Runs blastn for each sample against a local nucleotide BLAST database.
#   4) Writes per-sample BLAST results as TSV files.
#   5) Runs a Python script to summarise those TSVs into CSV outputs.
#
# Usage:
#   02_blast_cmv_summary.sh <FASTQ_DIR_OR_ARCHIVE> <BLAST_DB_DIR> <OUT_DIR> \
#       <THREADS>
#
# Inputs:
#   FASTQ_DIR_OR_ARCHIVE
#     - A directory containing .fastq/.fq (optionally gzipped), OR
#     - A supported archive containing FASTQs (.zip, .tar.*, etc.).
#
#   BLAST_DB_DIR
#     - Directory containing a nucleotide BLAST database named:
#       viral_nt_db.*
#
#   OUT_DIR
#     - Output directory where results will be written.
#
#   THREADS
#     - Number of BLAST threads to use (speeds up BLAST).
#
# Outputs:
#   OUT_DIR/fasta/
#     - FASTA files created from each input FASTQ.
#   OUT_DIR/blast_output/
#     - One TSV per sample: <sample>_blast.tsv
#   OUT_DIR/
#     - Summary CSVs written by summarize_blast_hits.py

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Helpers ----------------------------- #
log() {
    # Print a timestamped log message.
    printf '[%(%F %T)T] %s\n' -1 "$*"
}

need() {
    # Check a required command exists before running the pipeline.
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "❌ Error: required tool '$1' not found in PATH." >&2
        exit 127
    fi
}

TMP_EXTRACT=""

cleanup() {
    # Remove temporary extracted files if we created them.
    # This should never change the script exit code.
    if [[ -n "${TMP_EXTRACT:-}" && -d "${TMP_EXTRACT:-}" ]]; then
        rm -rf "$TMP_EXTRACT" || true
    fi
}

# Always run cleanup when the script exits (success or failure).
trap cleanup EXIT

# ------------------------------ Inputs ------------------------------ #
if [[ "$#" -ne 4 ]]; then
    echo "Usage: $0 <FASTQ_DIR_OR_ARCHIVE> <BLAST_DB_DIR> <OUT_DIR> <THREADS>" \
        >&2
    exit 1
fi

FASTQ_SRC="$1"      # Directory of FASTQs or an archive file containing FASTQs.
BLAST_DB_DIR="$2"   # Folder containing the BLAST database files.
OUT_DIR="$3"        # Where to write results.
THREADS="$4"        # Number of threads for BLAST.

# ----------------------- Check required tools ----------------------- #
need blastn
need seqtk
need zcat
need python3

# ------------------------ Create output folders --------------------- #
FASTA_DIR="${OUT_DIR}/fasta"
BLAST_OUT_DIR="${OUT_DIR}/blast_output"
mkdir -p "$FASTA_DIR" "$BLAST_OUT_DIR"

# -------------------------- Check BLAST DB -------------------------- #
# This script expects a nucleotide BLAST database called "viral_nt_db".
DB_PREFIX="${BLAST_DB_DIR%/}/viral_nt_db"

# BLAST databases are stored as multiple files with the same prefix.
# We check that at least one "viral_nt_db.*" file exists.
if ! ls "${DB_PREFIX}."* >/dev/null 2>&1; then
    echo "❌ Error: BLAST DB '${DB_PREFIX}' not found (expected files like " \
         "${DB_PREFIX}.n*)." >&2
    exit 1
fi

# -------------------- Handle FASTQ input type ----------------------- #
# The input FASTQs can be:
#   - a directory containing FASTQ files, OR
#   - an archive file (zip/tar/tar.gz/etc.) that contains FASTQ files.
WORK_FQ="$FASTQ_SRC"

if [[ -f "$FASTQ_SRC" ]]; then
    case "$FASTQ_SRC" in
        *.zip|*.tar|*.tar.gz|*.tgz|*.tar.bz2|*.tbz2|*.tar.xz|*.txz)
            log "Detected FASTQ archive: $FASTQ_SRC"
            TMP_EXTRACT="$(mktemp -d)"
            log "Extracting to $TMP_EXTRACT ..."

            # Use Python to unpack the archive into a temporary folder.
            SRC="$FASTQ_SRC" DST="$TMP_EXTRACT" python3 - <<'PY'
import os
import shutil

src = os.environ["SRC"]
dst = os.environ["DST"]
shutil.unpack_archive(src, dst)
PY

            WORK_FQ="$TMP_EXTRACT"
            ;;
        *)
            echo "❌ Error: FASTQ source is a file but not a supported archive." \
                >&2
            exit 1
            ;;
    esac
elif [[ -d "$FASTQ_SRC" ]]; then
    WORK_FQ="$FASTQ_SRC"
else
    echo "❌ Error: FASTQ source not found: $FASTQ_SRC" >&2
    exit 1
fi

# --------------------------- Run BLAST ------------------------------ #
# For each FASTQ file found:
#   1) Convert FASTQ -> FASTA
#   2) Run blastn and save a TSV output file
processed=0

while IFS= read -r -d '' fq; do
    [[ -e "$fq" ]] || continue

    # Create a base sample name from the FASTQ filename.
    fname="$(basename "$fq")"
    base="${fname%.gz}"
    base="${base%.*}"   # drop .fastq/.fq
    base="${base%.*}"   # extra safety if name contains multiple dots

    log "🔄 Processing: $base"

    # Convert FASTQ to FASTA (BLAST query input).
    fasta_file="${FASTA_DIR}/${base}.fasta"
    if [[ "$fq" == *.gz ]]; then
        zcat "$fq" | seqtk seq -a - > "$fasta_file"
    else
        seqtk seq -a "$fq" > "$fasta_file"
    fi

    # Run blastn against the database.
    log "🧬 Running BLAST for ${base}…"
    blastn \
        -query "$fasta_file" \
        -db "$DB_PREFIX" \
        -evalue 1e-5 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -num_threads "$THREADS" \
        > "${BLAST_OUT_DIR}/${base}_blast.tsv"

    processed=$((processed + 1))
done < <(
    find "$WORK_FQ" \
        -type f \
        -regextype posix-extended \
        -regex '.*\.(fastq|fq)(\.gz)?$' \
        -print0
)

# If we did not find any FASTQs, fail clearly.
if [[ "$processed" -eq 0 ]]; then
    echo "❌ Error: no FASTQ files found under '$WORK_FQ'." >&2
    exit 1
fi

log "✅ All BLAST runs complete. TSVs: ${BLAST_OUT_DIR}"
log "➡️  Summarizing BLAST hits into: ${OUT_DIR}"

# Run the summariser script on the TSV outputs.
python3 summarize_blast_hits.py "$BLAST_OUT_DIR" "$OUT_DIR"

log "✅ BLAST summary CSVs written to: ${OUT_DIR}"

# Exit success (cleanup still runs because of the trap).
exit 0
