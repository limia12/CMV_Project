#!/usr/bin/env bash
# build_snpeff_db.sh
#
# What this script does:
#   It prepares files so snpEff can add gene/impact annotations to a VCF.
#   It builds a small, local snpEff "database" called CMV_Custom using:
#     - a reference FASTA (the genome sequence)
#     - a GFF/GFF3 file (where genes/features are on that genome)
#
# Usage:
#   build_snpeff_db.sh <reference_fasta> <gff_or_gff3> <work_dir> \
#       [protein_fasta_or_blank]
#
# Notes:
#   - If the GFF is missing/empty, the script exits successfully (no build).
#   - If the GFF contains an embedded "##FASTA" block, it is ignored.
#   - If the GFF uses a different contig name than the FASTA, the script
#     tries to fix it (only when there is one contig).
#   - snpEff checks are disabled so small mismatches do not crash the build.

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Inputs ------------------------------ #
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <reference_fasta> <gff_or_gff3> <work_dir> " \
         "[protein_fasta_or_blank]" >&2
    exit 1
fi

ref_fa="$1"                # Reference genome FASTA.
gff_in="${2:-}"            # Gene annotation file (GFF/GFF3).
work="$3"                  # Folder to build the snpEff files in.
prot_fa="${4:-}"           # Optional protein FASTA (not required).
genome_id="CMV_Custom"     # Name of the snpEff database.

# If the GFF is missing or empty, skip the build and return success.
if [[ -z "$gff_in" || ! -s "$gff_in" ]]; then
    echo "[snpEff] No GFF provided. Skipping build."
    exit 0
fi

# Create folders snpEff expects.
mkdir -p "$work/data/$genome_id" "$work/logs"

# ------------------------ Read FASTA contig name -------------------- #
# Many tools need the "contig name" in the GFF to match the FASTA header.
# We read the first FASTA header and take its first word as the contig name.
if [[ ! -s "$ref_fa" ]]; then
    echo "ERROR: reference FASTA not found: $ref_fa" >&2
    exit 2
fi

ref_name="$(awk '/^>/{gsub(/^>/,"",$1); print $1; exit}' "$ref_fa)"
if [[ -z "$ref_name" ]]; then
    echo "ERROR: could not read contig name from FASTA header" >&2
    exit 2
fi
echo "[snpEff] Reference contig name: $ref_name"

# ---------------------- Read contig names in the GFF ---------------- #
# Column 1 in a GFF is the contig name.
# If the file contains a "##FASTA" section, stop before it.
mapfile -t gff_ids < <(awk '
    BEGIN{FS=OFS="\t"}
    /^##FASTA/ { exit }
    /^#/      { next }
    NF>=3 && $1!="" && $1!="." { print $1 }
' "$gff_in" | sort -u)

echo "[snpEff] GFF contig names: ${gff_ids[*]:-(none)}"

# This is where we will write a cleaned version of the GFF for snpEff.
clean_gff="$work/data/$genome_id/genes.gff"

# ------------------ Fix contig name mismatch (if simple) ------------ #
# If there is exactly one contig name in the GFF and it does not match the
# FASTA contig name, rewrite the GFF contig name to match the FASTA.
if (( ${#gff_ids[@]} == 1 )) && [[ "${gff_ids[0]}" != "$ref_name" ]]; then
    old="${gff_ids[0]}"
    echo "[snpEff] Renaming GFF contig '$old' -> '$ref_name'"

    awk -v tgt="$ref_name" -v old="$old" '
        BEGIN{FS=OFS="\t"}
        /^##FASTA/ { exit }

        # Update the "sequence-region" line if present.
        /^##sequence-region[ \t]/ {
            sub("^##sequence-region[ \t]+" old, "##sequence-region " tgt)
            print
            next
        }

        # Keep comment lines.
        /^#/ { print; next }

        # Rewrite contig name in feature lines.
        NF>=3 { $1=tgt; print }
    ' "$gff_in" > "$clean_gff"
else
    # Otherwise, just copy the GFF section and ignore any embedded FASTA.
    awk '
        /^##FASTA/ { exit }
        { print }
    ' "$gff_in" > "$clean_gff"
fi

# ------------------------ Copy files for snpEff --------------------- #
# snpEff expects these names inside data/<genome_id>/:
#   - sequences.fa  (the reference genome)
#   - genes.gff     (the annotation file)
cp -f "$ref_fa" "$work/data/$genome_id/sequences.fa"

# Optional: copy protein FASTA if provided (not required for this build).
if [[ -n "$prot_fa" && -s "$prot_fa" ]]; then
    cp -f "$prot_fa" "$work/data/$genome_id/protein.fa"
    echo "[snpEff] Protein FASTA copied (checks are disabled)."
elif [[ -n "$prot_fa" ]]; then
    echo "[snpEff] Protein FASTA path was provided but not found. Ignored."
fi

# -------------------------- Write snpEff config --------------------- #
# This local config tells snpEff where our "data" folder is.
config="$work/snpeff.config"
{
    echo "data.dir = $work/data"
    echo "$genome_id.genome : CMV_Custom"
} > "$config"

# Save useful values for other scripts (optional convenience).
echo "$genome_id" > "$work/GENOME_ID.txt"
echo "$config" > "$work/CONFIG_PATH.txt"

# ------------------------- Run snpEff command ----------------------- #
# This runs snpEff either from PATH, or from a jar file if SNPEFF_JAR is set.
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

# We turn off "-e" briefly so we can capture the real snpEff exit code
# even though we pipe output into tee (to save a log).
set +e
run_snpeff build -gff3 -v -c "$config" -noCheckCds -noCheckProtein \
    "$genome_id" | tee "$work/logs/build.log"
rc=${PIPESTATUS[0]}
set -e

# If snpEff failed, stop and point to the log.
if [[ $rc -ne 0 ]]; then
    echo "ERROR: snpEff build failed. See: $work/logs/build.log" >&2
    exit "$rc"
fi

echo "[snpEff] Build complete."
