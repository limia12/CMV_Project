#!/usr/bin/env bash
# alignment_bwa-mem.sh
# Usage:
#   alignment_bwa-mem.sh <fastq_dir_or_archive> <genome_index_dir_or_fasta> <output_root> [--threads N|N]
#
# Align paired-end FASTQs with BWA-MEM, sort/index BAMs, and write stats.

set -euo pipefail #prevents errors in a pipeline from being masked


# ---------------- helpers ----------------
# timestamped logging for debugging and monitoring
log() { printf '[%(%F %T)T] %s\n' -1 "$*"; }

# check for required tools
need() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Error: required tool '$1' not found in PATH." >&2
    exit 127
  fi
}

# If FASTQs are provided as an archive, extract to a temp folder and auto-clean it on exit
tmp_extract=""
cleanup() {
  if [[ -n "${tmp_extract:-}" && -d "${tmp_extract:-}" ]]; then
    rm -rf "$tmp_extract" || true
  fi
}
trap cleanup EXIT

# ---------------- args ----------------
# specifies required CLI arguments and prints usage if incorrect.
if [[ "$#" -lt 3 ]]; then
  echo "Usage: $0 <fastq_dir_or_archive> <genome_index_dir_or_fasta> <output_root> [--threads N|N]" >&2
  exit 1
fi

# defines CLI arguments and prints usage if incorrect.
input_path="$1"
genome_arg="$2"
output_root="$3"
threads=4

# support either "--threads N" or positional N
if [[ "${4:-}" == "--threads" ]]; then
  [[ "${5:-}" =~ ^[0-9]+$ ]] || { echo "Error: --threads requires an integer." >&2; exit 1; }
  threads="$5"
elif [[ "${4:-}" =~ ^[0-9]+$ ]]; then
  threads="$4"
fi

# Verify dependencies exist before starting alignment
# Stop immediately if any dependency is missing.
need bwa
need samtools
need python3

#Creates output directories for sorted BAMs and index/stats.
bam_dir="$output_root/BAM"
stats_dir="$bam_dir/BAM_Stats"
mkdir -p "$stats_dir"

# ---------------- normalize FASTQ input ----------------
#If input is already a folder, resolve to absolute path.
#If input is an archive, extract to a temp dir or reject unsupported formats.
if [[ -d "$input_path" ]]; then
  fastq_dir="$(realpath "$input_path")"
elif [[ -f "$input_path" ]]; then
  case "$input_path" in
    *.zip|*.tar|*.tar.gz|*.tgz|*.tar.bz2|*.tbz2|*.tar.xz|*.txz)
      log "Detected FASTQ archive: $input_path"
      tmp_extract="$(mktemp -d)"
      log "Extracting to $tmp_extract ..."
      python3 - <<'PY' "$input_path" "$tmp_extract"
import sys, shutil
src, dst = sys.argv[1], sys.argv[2]
shutil.unpack_archive(src, dst)
PY
      fastq_dir="$tmp_extract"
      ;;
    *)
      echo "Error: '$input_path' is a file but not a supported archive." >&2
      exit 1
      ;;
  esac
else
  echo "Error: '$input_path' not found." >&2
  exit 1
fi

# ---------------- find BWA index prefix ----------------
# Look for *.amb (present for bwa-mem2 and bwa index) and strip suffix
# If genome_arg is a FASTA: prefer an index with the same basename; otherwise search nearby for an index
# Stop early if we can't find BWA index files (expected *.amb/*.bwt from `bwa index`)
index_prefix=""
if [[ -d "$genome_arg" ]]; then
  if ! index_prefix="$(find "$genome_arg" -type f -name '*.amb' -print -quit)"; then
    index_prefix=""
  fi
  [[ -n "$index_prefix" ]] && index_prefix="${index_prefix%.amb}"
elif [[ -f "$genome_arg" ]]; then
  base="${genome_arg%.*}"
  if [[ -f "${base}.amb" || -f "${base}.bwt" ]]; then
    index_prefix="$base"
  else
    dir="$(dirname "$genome_arg")"
    if ! index_prefix="$(find "$dir" -type f -name '*.amb' -print -quit)"; then
      index_prefix=""
    fi
    [[ -n "$index_prefix" ]] && index_prefix="${index_prefix%.amb}"
  fi
fi

if [[ -z "$index_prefix" ]]; then
  echo "Error: could not locate BWA index in '$genome_arg'.
Expecting files like *.amb/*.bwt produced by 'bwa index'." >&2
  exit 1
fi

# ---------------- collect R1 files ----------------
# Find R1 FASTQs using common naming patterns and store them in an array (sorted)
# Abort if no R1 FASTQs were found (input naming may be unexpected)
cd "$fastq_dir"

mapfile -t R1S < <(
  find . -type f \( \
    -iname '*_R1.fastq' -o -iname '*_R1.fq' -o -iname '*_R1.fastq.gz' -o -iname '*_R1.fq.gz' -o \
    -iname '*_1.fastq'  -o -iname '*_1.fq'  -o -iname '*_1.fastq.gz'  -o -iname '*_1.fq.gz'  -o \
    -iname 'R1_*.fastq' -o -iname 'R1_*.fq' -o -iname 'R1_*.fastq.gz' -o -iname 'R1_*.fq.gz' \
  \) | sort
)

if [[ "${#R1S[@]}" -eq 0 ]]; then
  echo "Error: no FASTQs matched in '$fastq_dir' (expected *_R1.fastq[.gz], *_1.fastq[.gz], or R1_*.fastq[.gz])." >&2
  exit 1
fi

# ---------------- align each pair ----------------
# Loop over each detected R1 file; infer its mate (R2) and align as a paired-end sample
# Infer mate FASTQ (R2) filename from R1 using common naming conventions
# Skip if mate FASTQ is missing (warn but continue processing other samples)

for fq1 in "${R1S[@]}"; do
  fq1="${fq1#./}"

  if [[ "$fq1" =~ _R1\. ]]; then
    fq2="${fq1/_R1./_R2.}"
  elif [[ "$fq1" =~ _1\. ]]; then
    fq2="${fq1/_1./_2.}"
  elif [[ "$fq1" =~ R1_ ]]; then
    fq2="${fq1/R1_/R2_}"
  else
    fq2="${fq1/_R1/_R2}"
  fi

  if [[ ! -f "$fq2" ]]; then
    echo "⚠️  Paired file for $fq1 not found. Skipping." >&2
    continue
  fi

  # Derive sample ID from R1 filename by stripping extensions and the R1 marker (_R1 or _1)
  bn="$(basename "$fq1")"
  sample="$bn"
  sample="${sample%.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"
  sample="${sample%_R1}"
  sample="${sample%_1}"
  # (R1_ prefix style is uncommon for sample ID; left as-is)

  log "🔄 Aligning: $sample (threads=$threads)"
  sorted_bam="$bam_dir/${sample}_sorted.bam"

# Align paired-end reads with BWA-MEM, then stream directly into samtools sort to produce a coordinate-sorted BAM
# -M marks shorter split hits as secondary (compatibility with some downstream pipelines)
  bwa mem -M -t "$threads" "$index_prefix" "$fq1" "$fq2" \
    | samtools sort -@ "$threads" -o "$sorted_bam" -
    
# Index BAM for random access, then write alignment QC summaries:
# flagstat = high-level mapping summary; stats = detailed alignment statistics
  samtools index -@ "$threads" "$sorted_bam"
  samtools flagstat -@ "$threads" "$sorted_bam" > "$stats_dir/${sample}_flagstat.txt"
  samtools stats    -@ "$threads" "$sorted_bam" > "$stats_dir/${sample}_stats.txt"

  log "✅ Finished: $sample"
done

log "🎉 All alignments complete. Sorted BAMs in '$bam_dir'."
exit 0
