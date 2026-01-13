#!/usr/bin/env bash
# filtering.sh
# Usage: filtering.sh <input_dir_or_archive> <output_dir> <quality_threshold> <num_cores>
# Trims paired-end FASTQs with Trim Galore and runs FastQC.
# Accepts a directory of FASTQs or a supported archive (zip/tar*).

set -euo pipefail #prevents errors in a pipeline from being masked

# ---------- helpers ----------
# timestamped logging for debugging and monitoring
log() { printf '[%(%F %T)T] %s\n' -1 "$*"; }

# check for required tools
need() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Error: required tool '$1' not found in PATH." >&2
    exit 127
  fi
}

# If input is an archive, extract to a temp dir
# Cleanup function to remove temp dir on exit.
tmpdir=""
cleanup() {
  # Never let cleanup alter exit status.
  if [[ -n "${tmpdir:-}" && -d "${tmpdir:-}" ]]; then
    rm -rf "$tmpdir" || true
  fi
}
trap cleanup EXIT

# ---------- args ----------
# defines CLI arguments and prints usage if incorrect.
if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 <input_dir_or_archive> <output_dir> <quality_threshold> <num_cores>" >&2
  exit 1
fi

input=$1
output_dir=$2
quality_threshold=$3
num_cores=$4

# Validation - Ensures parameters are integers and cores are more than 0.
[[ "$quality_threshold" =~ ^[0-9]+$ ]] || { echo "Error: quality_threshold must be an integer." >&2; exit 1; }
[[ "$num_cores" =~ ^[0-9]+$ ]]        || { echo "Error: num_cores must be an integer." >&2; exit 1; }
[[ "$num_cores" -ge 1 ]]              || { echo "Error: num_cores must be >= 1." >&2; exit 1; }

# Ensures required rools are available stops immediately if any dependency is missing. 
need trim_galore
need fastqc
need python3

# If input is an archive, unpack to a temp folder so downstream steps always work on a directory
workdir="$input"
if [[ -f "$input" ]]; then
  case "$input" in
    *.zip|*.tar|*.tar.gz|*.tgz|*.tar.bz2|*.tbz2|*.tar.xz|*.txz)
      tmpdir="$(mktemp -d)"
      log "Extracting archive: $input -> $tmpdir"
      python3 - <<'PY' "$input" "$tmpdir"
import sys, shutil
shutil.unpack_archive(sys.argv[1], sys.argv[2])
PY
      workdir="$tmpdir"
      ;;
    *)
      echo "Error: '$input' is a file but not a supported archive." >&2
      exit 1
      ;;
  esac
fi

# Stops if the working directory does not exist
if [[ ! -d "$workdir" ]]; then
  echo "Error: Input directory '$workdir' does not exist." >&2
  exit 1
fi

# Creates output directories and removes previous run outputs if they exist
filtered_dir="$output_dir/filtered_reads"
fastqc_dir="$output_dir/fastqc_results"

log "Setting up output directories…"
mkdir -p "$filtered_dir" "$fastqc_dir"
rm -rf "$filtered_dir"/* "$fastqc_dir"/* 2>/dev/null || true

log "Starting paired-end trimming using Trim Galore…"
log "Quality threshold: $quality_threshold | Cores: $num_cores"

#Export vars so each parallel xargs worker can access output directory + trimming parameters

export filtered_dir quality_threshold

# ---------- trim in parallel (by R1) ----------
# Find all R1 files and process each pair in parallel up to num_cores samples at once
find "$workdir" -type f \( \
   -iname '*_R1.fastq' -o -iname '*_R1.fastq.gz' -o -iname '*_R1.fq' -o -iname '*_R1.fq.gz' -o \
   -iname '*_1.fastq'  -o -iname '*_1.fastq.gz'  -o -iname '*_1.fq'  -o -iname '*_1.fq.gz' \
\) -print0 \
| xargs -0 -n1 -P "$num_cores" bash -c '
  set -euo pipefail
  r1="$1"

  # Find the matching R2 name from R1 filename 
  if [[ "$r1" =~ _R1\. ]]; then
    r2="${r1/_R1./_R2.}"
  elif [[ "$r1" =~ _1\. ]]; then
    r2="${r1/_1./_2.}"
  elif [[ "$r1" =~ R1_ ]]; then
    r2="${r1/R1_/R2_}"
  else
    r2="${r1/_R1/_R2}"
  fi
  
  # Extract sample ID from R1 filename skips unpaired files rather than failing whole batch
  base="$(basename "$r1")"
  sample="$(echo "$base" | sed -E "s/(_R1|_1)\\.(fastq|fq)(\\.gz)?$//I")"

  if [[ ! -f "$r2" ]]; then
    echo "⚠️  Paired file for $r1 not found. Skipping $sample." >&2
    exit 0
  fi

  echo "Processing $sample…"

  # Run Trim Galore with adapter and quality trimming with Trim Galore
  if ! trim_galore --quality "$quality_threshold" --paired --cores 1 --output_dir "$filtered_dir" "$r1" "$r2"; then
    echo "❌ trim_galore failed for $sample" >&2
    exit 1
  fi

  # Track whether outputs are gzipped so rename logic matches the correct filenames
  gz=""
  [[ "$r1" == *.gz ]] && gz=".gz"

  # standardize output filenames so downstream steps can find them easily.
  if [[ -f "$filtered_dir/${sample}_R1_val_1.fq$gz" ]]; then
    mv -f "$filtered_dir/${sample}_R1_val_1.fq$gz" "$filtered_dir/${sample}_R1.fastq$gz"
  elif [[ -f "$filtered_dir/${sample}_R1_val_1.fq" ]]; then
    mv -f "$filtered_dir/${sample}_R1_val_1.fq" "$filtered_dir/${sample}_R1.fastq"
  fi

  if [[ -f "$filtered_dir/${sample}_R2_val_2.fq$gz" ]]; then
    mv -f "$filtered_dir/${sample}_R2_val_2.fq$gz" "$filtered_dir/${sample}_R2.fastq$gz"
  elif [[ -f "$filtered_dir/${sample}_R2_val_2.fq" ]]; then
    mv -f "$filtered_dir/${sample}_R2_val_2.fq" "$filtered_dir/${sample}_R2.fastq"
  fi
' _

# confirm outputs exist
# Sanity check - ensure trimming actually produced FASTQs (prevents FastQC running on empty folder)
if ! find "$filtered_dir" -type f \
        \( -iname '*.fastq' -o -iname '*.fq' -o -iname '*.fastq.gz' -o -iname '*.fq.gz' \) \
        -print -quit | grep -q . ; then
  echo "Error: no trimmed FASTQs were produced. Check input naming and trim_galore logs." >&2
  exit 1
fi

log "Trimming complete. Starting FastQC…"

# Run FastQC on trimmed reads. Parallel across FASTQ files (one thread per FastQC job).
find "$filtered_dir" -type f \( -iname '*.fastq' -o -iname '*.fq' -o -iname '*.fastq.gz' -o -iname '*.fq.gz' \) -print0 \
| xargs -0 -n1 -P "$num_cores" fastqc -t 1 --outdir "$fastqc_dir"

log "FastQC complete."
log "Filtered reads saved in: $filtered_dir"
log "FastQC results saved in: $fastqc_dir"

exit 0
