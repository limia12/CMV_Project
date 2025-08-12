#!/bin/bash
# variant_calling_vaf.sh
# Usage: variant_calling_vaf.sh <bam_dir> <reference_fasta> <output_dir> [threads]
# - Calls variants with LoFreq and produces a VAF table using INFO fields (DP, AF, DP4).

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <bam_dir> <reference_fasta> <output_dir> [threads]" >&2
    exit 1
fi

bam_dir="$1"
ref_fa="$2"
out_dir="$3"
threads="${4:-4}"

# Remove old outputs if they exist
rm -rf "$out_dir/vcf" "$out_dir/vaf_tables"
rm -f  "$out_dir"/*_indelqual.bam "$out_dir"/*_indelqual.bam.bai

# Recreate directories
mkdir -p "$out_dir/vcf" "$out_dir/vaf_tables"

# Index reference if needed
if [ ! -f "${ref_fa}.fai" ]; then
    samtools faidx "$ref_fa"
fi

for bam in "$bam_dir"/*sorted.bam; do
    [ -f "$bam" ] || continue
    sample="$(basename "$bam" _sorted.bam)"
    echo "[INFO] Processing $sample..."

    # Ensure BAM is indexed
    if [ ! -f "${bam}.bai" ]; then
        samtools index -@ "$threads" "$bam"
    fi

    # Add indel qualities required by LoFreq
    indel_bam="$out_dir/${sample}_indelqual.bam"
    lofreq indelqual --dindel -f "$ref_fa" -o "$indel_bam" "$bam"
    samtools index "$indel_bam"

    # Call variants
    vcf_out="$out_dir/vcf/${sample}.vcf"
    lofreq call-parallel --pp-threads "$threads" -f "$ref_fa" -o "$vcf_out" "$indel_bam"

    # Build VAF table from INFO tags
    vaf_out="$out_dir/vaf_tables/${sample}_vaf.tsv"
    echo -e "CHROM\tPOS\tREF\tALT\tDP\tAD_ALT\tAF" > "$vaf_out"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t%INFO/DP4\n' "$vcf_out" \
    | awk 'BEGIN{OFS="\t"}{
        dp = ($5 == "." || $5 == "" ? 0 : $5)
        alt=0
        if ($7 != "." && $7 != "") {
            split($7, d, ",");
            if (length(d) >= 4) { alt = d[3] + d[4] }
        }
        af = ($6 == "." || $6 == "" ? (dp>0 ? alt/dp : 0) : $6)
        print $1,$2,$3,$4,dp,alt,af
    }' >> "$vaf_out"

done

echo "[INFO] Variant calling + VAF extraction complete."
