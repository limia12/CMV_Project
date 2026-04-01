#!/usr/bin/env bash
# build_snpeff_db.sh
#
# Build a local snpEff database from a reference FASTA + GFF3 for CMV.
#
# Usage:
#   bash build_snpeff_db.sh <reference_fasta> <gff_or_gff3> <work_dir> [protein_fasta_or_blank]
#
# Notes:
# - If the GFF is missing/empty, the script exits 0 (skip build).
# - If the GFF contains a "##FASTA" block, it is ignored.
# - If there is exactly ONE contig in the GFF and it doesn't match the FASTA contig,
#   the GFF contig is rewritten to match the FASTA contig.
# - snpEff checks are disabled to tolerate small annotation quirks.

set -euo pipefail

# ------------------------------ Inputs ------------------------------ #
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <reference_fasta> <gff_or_gff3> <work_dir> [protein_fasta_or_blank]" >&2
  exit 1
fi

ref_fa="$1"
gff_in="$2"
work="$3"
prot_fa="${4:-}"
genome_id="CMV_Custom"

# Skip build if no GFF
if [[ -z "${gff_in:-}" || ! -s "$gff_in" ]]; then
  echo "[snpEff] No GFF provided (or empty). Skipping build."
  exit 0
fi

if [[ ! -s "$ref_fa" ]]; then
  echo "ERROR: reference FASTA not found or empty: $ref_fa" >&2
  exit 2
fi

mkdir -p "$work/data/$genome_id" "$work/logs"

# ------------------------ Read FASTA contig name -------------------- #
ref_name="$(
  awk '
    /^>/ {
      gsub(/^>/,"",$1);
      print $1;
      exit
    }
  ' "$ref_fa"
)"

if [[ -z "$ref_name" ]]; then
  echo "ERROR: could not read contig name from FASTA header: $ref_fa" >&2
  exit 2
fi
echo "[snpEff] Reference contig name: $ref_name"

# ---------------------- Read contig names in the GFF ---------------- #
# Stop before embedded FASTA, ignore comments, collect col1
gff_ids_tmp="$work/logs/gff_contigs.txt"
awk '
  BEGIN{FS=OFS="\t"}
  /^##FASTA/ { exit }
  /^#/       { next }
  NF>=3 && $1!="" && $1!="." { print $1 }
' "$gff_in" | sort -u > "$gff_ids_tmp"

# Read into bash array safely
mapfile -t gff_ids < "$gff_ids_tmp"
if (( ${#gff_ids[@]} == 0 )); then
  echo "[snpEff] WARNING: No contig IDs detected in GFF (feature lines). Build may fail."
else
  echo "[snpEff] GFF contig names: ${gff_ids[*]}"
fi

clean_gff="$work/data/$genome_id/genes.gff"

# ------------------ Fix contig name mismatch (if simple) ------------ #
if (( ${#gff_ids[@]} == 1 )) && [[ "${gff_ids[0]}" != "$ref_name" ]]; then
  old="${gff_ids[0]}"
  echo "[snpEff] Renaming GFF contig '$old' -> '$ref_name'"

  awk -v tgt="$ref_name" -v old="$old" '
    BEGIN{FS=OFS="\t"}
    /^##FASTA/ { exit }

    /^##sequence-region[ \t]/ {
      sub("^##sequence-region[ \t]+" old, "##sequence-region " tgt)
      print
      next
    }

    /^#/ { print; next }

    NF>=3 { $1=tgt; print; next }

    { print }
  ' "$gff_in" > "$clean_gff"
else
  awk '
    /^##FASTA/ { exit }
    { print }
  ' "$gff_in" > "$clean_gff"
fi

# ------------------------ Copy files for snpEff --------------------- #
cp -f "$ref_fa" "$work/data/$genome_id/sequences.fa"

if [[ -n "${prot_fa:-}" && -s "$prot_fa" ]]; then
  cp -f "$prot_fa" "$work/data/$genome_id/protein.fa"
  echo "[snpEff] Protein FASTA copied."
elif [[ -n "${prot_fa:-}" ]]; then
  echo "[snpEff] Protein FASTA path provided but not found/empty. Ignored."
fi

# -------------------------- Write snpEff config --------------------- #
config="$work/snpeff.config"
{
  echo "data.dir = $work/data"
  echo "$genome_id.genome : CMV_Custom"
} > "$config"

echo "$genome_id" > "$work/GENOME_ID.txt"
echo "$config" > "$work/CONFIG_PATH.txt"

# ------------------------- Run snpEff command ----------------------- #
run_snpeff() {
  if command -v snpEff >/dev/null 2>&1; then
    snpEff "$@"
  elif [[ -n "${SNPEFF_JAR:-}" && -f "${SNPEFF_JAR:-}" ]]; then
    java -Xmx4g -jar "$SNPEFF_JAR" "$@"
  else
    echo "ERROR: snpEff not found. Install it or set SNPEFF_JAR." >&2
    exit 3
  fi
}

echo "[snpEff] Building database: $genome_id"
echo "[snpEff] Using:"
echo "         FASTA: $work/data/$genome_id/sequences.fa"
echo "         GFF  : $work/data/$genome_id/genes.gff"
echo "         CONF : $config"

set +e
run_snpeff build -gff3 -v -c "$config" -noCheckCds -noCheckProtein "$genome_id" \
  2>&1 | tee "$work/logs/build.log"
rc=${PIPESTATUS[0]}
set -e

if [[ $rc -ne 0 ]]; then
  echo "ERROR: snpEff build failed (exit $rc). See log: $work/logs/build.log" >&2
  exit "$rc"
fi

echo "[snpEff] Build complete."
