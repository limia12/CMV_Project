#!/usr/bin/env bash
# annotate_vcfs_snpeff.sh
#
# What this script does:
#   It takes VCF files in a folder and runs snpEff to add annotations.
#   "Annotation" here means snpEff adds extra INFO fields describing what a
#   variant might do in/near genes (based on the snpEff database you built).
#
# Usage:
#   annotate_vcfs_snpeff.sh <vcf_dir> <snpeff_work_dir> [threads]
#
# Inputs:
#   vcf_dir         Folder containing VCF files (*.vcf or *.vcf.gz).
#   snpeff_work_dir Folder that contains snpEff build outputs and config.
#                  This is the same folder used by build_snpeff_db.sh.
#   threads         Optional threads for bgzip compression (default: 4).
#
# Outputs:
#   For each input VCF:
#     - If input is .vcf.gz  -> writes .annot.vcf.gz + creates a .tbi index
#     - If input is .vcf     -> writes .annot.vcf
#
# Notes:
#   - If no VCFs are found, the script exits successfully.
#   - This script supports snpEff installed on PATH or via SNPEFF_JAR.

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Inputs ------------------------------ #
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <vcf_dir> <snpeff_work_dir> [threads]" >&2
    exit 1
fi

vcf_dir="$1"          # Directory containing VCF files to annotate.
work="$2"             # snpEff work directory (where config + genome ID live).
threads="${3:-4}"     # Threads used for bgzip compression if available.

# ----------------------- Read genome + config ----------------------- #
# These files are written by the snpEff build script. If they are missing,
# we fall back to sensible defaults.
genome_id="$(tr -d '\r' < "$work/GENOME_ID.txt" 2>/dev/null || echo CMV_Custom)"
config_path="$(tr -d '\r' < "$work/CONFIG_PATH.txt" 2>/dev/null || echo "$work/snpeff.config")"

# Confirm the config file exists before continuing.
if [[ ! -f "$config_path" ]]; then
    echo "ERROR: snpEff config not found: $config_path" >&2
    exit 2
fi

echo "[snpEff] Using genome: $genome_id"
echo "[snpEff] Using config: $config_path"

# ------------------------ Run snpEff safely ------------------------- #
# This helper runs snpEff either:
#   - from the snpEff command (if installed), OR
#   - from a jar file if SNPEFF_JAR is set.
run_snpeff() {
    if command -v snpEff >/dev/null 2>&1; then
        snpEff "$@"
    elif [[ -n "${SNPEFF_JAR:-}" && -f "${SNPEFF_JAR:-}" ]]; then
        java -Xmx4g -jar "$SNPEFF_JAR" "$@"
    else
        echo "ERROR: snpEff not found. Install snpEff or set SNPEFF_JAR." >&2
        exit 3
    fi
}

# -------------------- Set up compression + indexing ----------------- #
# If bgzip exists, we use it to create compressed VCFs and then index them
# with tabix. This makes large annotated VCFs faster to open later.
bgzip_cmd=(bgzip -c)
tabix_cmd=(tabix -f -p vcf)

# If bgzip supports threads, use them.
if command -v bgzip >/dev/null 2>&1; then
    bgzip_cmd=(bgzip -@ "$threads" -c)
fi

# --------------------------- Find input VCFs ------------------------ #
# nullglob avoids literal patterns if no files match.
shopt -s nullglob

# Accept either uncompressed or gzipped VCF files.
vfs=( "$vcf_dir"/*.vcf "$vcf_dir"/*.vcf.gz )

# If there are no VCFs, exit successfully.
if [[ ${#vfs[@]} -eq 0 ]]; then
    echo "[snpEff] No VCFs found in $vcf_dir"
    exit 0
fi

# --------------------------- Annotate VCFs -------------------------- #
for v in "${vfs[@]}"; do
    base="$(basename "$v")"

    # If input is gzipped, keep output gzipped and create a tabix index.
    if [[ "$v" == *.gz ]]; then
        out="$vcf_dir/${base%.vcf.gz}.annot.vcf.gz"
        echo "[snpEff] Annotating $base -> $(basename "$out")"

        # snpEff writes to stdout. We compress that output into .vcf.gz.
        run_snpeff ann -c "$config_path" -v "$genome_id" -noLog "$v" \
            | "${bgzip_cmd[@]}" > "$out"

        # Create a .tbi index for the gzipped VCF.
        "${tabix_cmd[@]}" "$out"

    # If input is not gzipped, write a plain annotated VCF.
    else
        out="$vcf_dir/${base%.vcf}.annot.vcf"
        echo "[snpEff] Annotating $base -> $(basename "$out")"

        run_snpeff ann -c "$config_path" -v "$genome_id" -noLog "$v" > "$out"
    fi
done

echo "[snpEff] Annotation complete."
