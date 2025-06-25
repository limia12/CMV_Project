#!/bin/bash

# Usage check
if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <forward.fastq> <reverse.fastq> <blast_db_directory> <num_threads> <output_dir>"
  exit 1
fi

# Input arguments
FORWARD_FASTQ="$1"
REVERSE_FASTQ="$2"
BLAST_DB_DIR="$3"
NUM_THREADS="$4"
OUTPUT_DIR="$5"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Step 1: Identify BLAST database prefixes
echo "Searching for BLAST databases in: $BLAST_DB_DIR"
db_files=($(find "$BLAST_DB_DIR" -type f -name '*_db*'))

declare -A db_prefixes
for filepath in "${db_files[@]}"; do
  filename=$(basename "$filepath")
  prefix="${filename%%.*}"
  db_prefixes["$prefix"]=1
done

if [ "${#db_prefixes[@]}" -eq 0 ]; then
  echo "No BLAST databases found containing '_db' in $BLAST_DB_DIR"
  exit 1
fi

echo "BLAST database prefixes found:"
for dbp in "${!db_prefixes[@]}"; do
  echo "  $dbp"
done

# Step 2: Convert FASTQ to FASTA
echo "Converting FASTQ files to FASTA..."
seqtk seq -a "$FORWARD_FASTQ" > "$OUTPUT_DIR/reads_forward.fasta"
seqtk seq -a "$REVERSE_FASTQ" > "$OUTPUT_DIR/reads_reverse.fasta"

# Step 3: Run BLAST and collect results
COMBINED_OUT="$OUTPUT_DIR/blast_combined_results.tsv"
> "$COMBINED_OUT"  # Clear output file

SAMPLE_FWD=$(basename "$FORWARD_FASTQ")
SAMPLE_REV=$(basename "$REVERSE_FASTQ")

for db_prefix in "${!db_prefixes[@]}"; do
  DB_PATH="$BLAST_DB_DIR/$db_prefix"

  echo "Running blastn for $SAMPLE_FWD against $db_prefix..."
  blastn -query "$OUTPUT_DIR/reads_forward.fasta" \
         -db "$DB_PATH" \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
         -evalue 1e-5 \
         -num_threads "$NUM_THREADS" | \
  awk -v sample="$SAMPLE_FWD" '{print sample "\t" $0}' >> "$COMBINED_OUT"

  echo "Running blastn for $SAMPLE_REV against $db_prefix..."
  blastn -query "$OUTPUT_DIR/reads_reverse.fasta" \
         -db "$DB_PATH" \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
         -evalue 1e-5 \
         -num_threads "$NUM_THREADS" | \
  awk -v sample="$SAMPLE_REV" '{print sample "\t" $0}' >> "$COMBINED_OUT"
done

echo "BLAST runs complete. Results saved to $COMBINED_OUT"

# Step 4: Summarize results
SUMMARY_FILE="$OUTPUT_DIR/blast_summary.tsv"
TEMP_SUMMARY="$OUTPUT_DIR/blast_summary_unsorted.tsv"

# Write header
echo -e "Sample\tAccession\tHits\tAvg_Identity(%)\tVirus_Name" > "$SUMMARY_FILE"

# Parse results and compute statistics
awk '
{
  sample = $1            # prepended sample name
  acc = $3               # sseqid (accession)
  identity = $4          # percent identity
  virus_name = ""

  # Reconstruct virus name (stitle)
  for(i=14; i<=NF; i++) virus_name = virus_name $i (i==NF ? "" : " ")
  sub(acc "[[:space:]]+", "", virus_name)

  hits[sample, acc]++
  identity_sum[sample, acc] += identity
  virus_names[sample, acc] = virus_name
}
END {
  for (key in hits) {
    split(key, arr, SUBSEP)
    sample = arr[1]
    acc = arr[2]
    avg_identity = identity_sum[key] / hits[key]
    printf "%s\t%s\t%d\t%.2f\t%s\n", sample, acc, hits[key], avg_identity, virus_names[key]
  }
}
' "$COMBINED_OUT" > "$TEMP_SUMMARY"

# Sort by number of hits descending and finalize summary
sort -k3,3nr "$TEMP_SUMMARY" >> "$SUMMARY_FILE"
rm "$TEMP_SUMMARY"

echo "Summary written to $SUMMARY_FILE"
