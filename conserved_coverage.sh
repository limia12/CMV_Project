#!/bin/bash

# conserved_coverage.sh
# Usage: ./conserved_coverage.sh -i <input_bam_dir> -r <reference_fasta> -o <output_csv>

usage() {
  echo "Usage: $0 -i <input_bam_dir> -r <reference_fasta> -o <output_csv>"
  exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -i|--input) BAM_DIR="$2"; shift ;;
    -r|--reference) REFERENCE="$2"; shift ;;
    -o|--output) OUTPUT="$2"; shift ;;
    *) usage ;;
  esac
  shift
done

if [[ -z "$BAM_DIR" || -z "$REFERENCE" || -z "$OUTPUT" ]]; then
  usage
fi

# Check required tools
for tool in samtools awk; do
  command -v $tool >/dev/null 2>&1 || { echo "$tool is required but not installed." >&2; exit 1; }
done

# Index reference
if [[ ! -f "$REFERENCE.fai" ]]; then
  echo "Indexing reference FASTA..."
  samtools faidx "$REFERENCE"
fi

# Output header with new Depth column
echo "Sample,Gene,Percent_Alignment,Average_Depth" > "$OUTPUT"

regions=(
  "IE1:merlinGenome:172329-174090"
  "UL54:merlinGenome:78194-81922"
  "gpUL55:merlinGenome:82066-84789"
  "UL75:merlinGenome:109224-111452"
  "UL97:merlinGenome:141798-143921"
  "UL83:merlinGenome:120656-122341"
)

for bam in "$BAM_DIR"/*.bam; do
  [[ -f "$bam" ]] || continue
  sample=$(basename "$bam" .bam)

  for region in "${regions[@]}"; do
    IFS=":" read -r gene chrom coords <<< "$region"
    IFS="-" read -r start end <<< "$coords"
    region_length=$((end - start + 1))

    # Count covered bases (depth > 0)
    covered_bases=$(samtools depth -r "$chrom:$start-$end" "$bam" | awk '$3 > 0 {c++} END {print c+0}')

    # Calculate percent covered
    percent=$(awk -v cov="$covered_bases" -v len="$region_length" 'BEGIN { printf "%.2f", (cov/len)*100 }')

    # Calculate average depth over the region
    # If no coverage lines, avoid division by zero by setting depth=0
    avg_depth=$(samtools depth -r "$chrom:$start-$end" "$bam" | awk '{sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print 0}')

    echo "$sample,$gene,$percent,$avg_depth" >> "$OUTPUT"
  done
done

echo "Conserved gene coverage report generated at: $OUTPUT"
