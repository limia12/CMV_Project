#!/bin/bash
# variant_calling_vaf.sh
#
# Purpose:
#   Call variants with LoFreq for each "*sorted.bam" file in a directory.
#   Then write a simple TSV per sample with depth (DP), ALT depth (AD_ALT),
#   and allele frequency (AF / VAF).
#
# Usage:
#   variant_calling_vaf.sh <bam_dir> <ref_fasta(.gz)> <output_dir> [threads]
#
# Output:
#   <output_dir>/vcf/<sample>.vcf
#   <output_dir>/vaf_tables/<sample>_vaf.tsv

set -euo pipefail
# Exit on error (-e), unset variables (-u), or failed pipeline (pipefail).

# ------------------------------ Inputs ------------------------------ #
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <bam_dir> <reference_fasta_or_fasta.gz> <output_dir> [threads]" \
        >&2
    exit 1
fi

bam_dir="$1"        # Directory containing input BAMs named "*sorted.bam".
ref_in="$2"         # Reference FASTA (plain or gzipped).
out_dir="$3"        # Base output directory.
threads="${4:-4}"   # Threads for tools that support it (default: 4).

# -------------------------- Output folders -------------------------- #
# Remove old outputs so results are not mixed across runs.
rm -rf "$out_dir/vcf" "$out_dir/vaf_tables"
rm -f "$out_dir"/*_indelqual.bam "$out_dir"/*_indelqual.bam.bai

# Create fresh output directories.
mkdir -p "$out_dir/vcf" "$out_dir/vaf_tables"

# --------------------------- Reference FASTA ------------------------ #
# LoFreq typically works best with an uncompressed FASTA. If the reference
# is gzipped, decompress it to a temporary folder for this run.
ref_fa="$ref_in"
tmp_ref=""

if [[ "$ref_in" == *.gz ]]; then
    tmp_ref="$(mktemp -d)"
    ref_fa="$tmp_ref/$(basename "${ref_in%.gz}")"
    echo "[INFO] Decompressing reference: $ref_in -> $ref_fa"
    gzip -cd "$ref_in" > "$ref_fa"
fi

# Ensure the FASTA index exists (required by many downstream tools).
if [ ! -f "${ref_fa}.fai" ]; then
    samtools faidx "$ref_fa"
fi

# If no BAM files match the glob, do not treat the pattern as a literal.
shopt -s nullglob

# ------------------------- Per-sample loop -------------------------- #
for bam in "$bam_dir"/*sorted.bam; do
    [ -f "$bam" ] || continue

    # Use the BAM name to define the sample ID.
    sample="$(basename "$bam" _sorted.bam)"
    echo "[INFO] Processing sample: $sample"

    # Index the input BAM if needed (creates a .bai file).
    if [ ! -f "${bam}.bai" ]; then
        samtools index -@ "$threads" "$bam"
    fi

    # Add indel-quality scores for improved LoFreq calling accuracy.
    indel_bam="$out_dir/${sample}_indelqual.bam"
    lofreq indelqual --dindel -f "$ref_fa" -o "$indel_bam" "$bam"
    samtools index "$indel_bam"

    # Call variants with LoFreq and write a per-sample VCF.
    vcf_out="$out_dir/vcf/${sample}.vcf"
    lofreq call-parallel --pp-threads "$threads" -f "$ref_fa" -o "$vcf_out" \
        "$indel_bam"

    # Create a TSV header for VAF output.
    vaf_out="$out_dir/vaf_tables/${sample}_vaf.tsv"
    echo -e "CHROM\tPOS\tREF\tALT\tDP\tAD_ALT\tAF" > "$vaf_out"

    # Extract fields from INFO:
    #   DP  = total depth at site
    #   AF  = allele frequency (if present)
    #   DP4 = ref_fwd, ref_rev, alt_fwd, alt_rev
    #
    # AD_ALT is computed as alt_fwd + alt_rev (DP4 fields 3 and 4).
    # If AF is missing, AF is computed as AD_ALT / DP (if DP > 0).
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t%INFO/DP4\n' \
        "$vcf_out" \
    | awk 'BEGIN{OFS="\t"}{
        dp = ($5 == "." || $5 == "" ? 0 : $5)

        alt = 0
        if ($7 != "." && $7 != "") {
            split($7, d, ",")
            if (length(d) >= 4) {
                alt = d[3] + d[4]
            }
        }

        af = ($6 == "." || $6 == "" ? (dp > 0 ? alt / dp : 0) : $6)
        print $1, $2, $3, $4, dp, alt, af
    }' >> "$vaf_out"
done

# ------------------------------ Cleanup ----------------------------- #
# Remove the temporary reference folder if one was created.
if [ -n "$tmp_ref" ]; then
    rm -rf "$tmp_ref"
fi

echo "[INFO] Variant calling and VAF table generation complete."
