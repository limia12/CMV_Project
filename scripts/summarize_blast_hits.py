#!/usr/bin/env python3
"""
summarize_blast_hits.py

This script reads BLAST tabular outputs (TSV files) and writes simple summary
CSV files that are easier to inspect and plot.

It is designed to work with TSV files created using this BLAST outfmt:

    outfmt 6:
      qseqid sseqid pident length mismatch gapopen qstart qend
      sstart send evalue bitscore stitle

Expected input filenames
------------------------
The script expects BLAST TSV files named like:

    <something>_blast.tsv

This matches the output produced by the pipeline script that runs BLAST.

What the script produces
------------------------
For each sample base name, the script can write:

1) <sample>_blast_top5.csv
   - Shows the top 5 species (by hit count, then identity) for that sample.

2) <sample>_blast_heatmap.csv
   - Only for CMV-like hits (matching keywords such as "Human herpesvirus 5").
   - Contains intervals (start/end) and identity, suitable for plotting a
     simple "where do reads hit the CMV genome" heatmap.

3) blast_heatmap.csv
   - A combined file containing heatmap intervals for all samples.

Important notes 
----------------
- A BLAST "hit" is one row in the TSV output.
- pident means "percent identity" (how similar the alignment is).
- "Species" is taken from the first part of the BLAST subject title (stitle)
  before the first comma. 

"""

import argparse
import glob
import os
import re
from collections import defaultdict

import pandas as pd


def clean_species_name(stitle):
    """
    Convert a BLAST subject title into a simple species label.

    Parameters
    ----------
    stitle : str
        The subject title field from BLAST (stitle).

    Returns
    -------
    str
        The part of stitle before the first comma, stripped of whitespace.

    Why this exists
    ---------------
    BLAST titles can be long (strain, isolate, gene info, etc.).
    For quick summaries, we want a short "species-like" label.
    """
    return stitle.split(",")[0].strip()


def parse_sample_base_and_read(filename):
    """
    Extract a sample base name and (optional) read label from a filename.

    Parameters
    ----------
    filename : str
        Path to a BLAST TSV file. The script uses the basename only.

    Returns
    -------
    tuple[str, str]
        (sample_base, read)

        sample_base:
          A simplified sample ID.

        read:
          One of "R1", "R2", or "unknown".

    How the parsing works
    ---------------------
    1) Remove the folder path.
    2) Remove the "_blast" suffix and everything after it.
    3) Detect if the remaining name ends with "_R1" or "_R2".
       - If it ends cleanly, split into base + read.
       - Otherwise, try to remove "_R1" or "_R2" if found somewhere inside.
       - If no read tag is found, return "unknown".
    """
    basename = os.path.basename(filename)
    sample_part = basename.split("_blast")[0]

    if re.search(r"_R[12]$", sample_part):
        return sample_part.rsplit("_", 1)[0], sample_part.rsplit("_", 1)[1]
    if "_R1" in sample_part or "_R2" in sample_part:
        return sample_part.replace("_R1", "").replace("_R2", ""), "unknown"
    return sample_part, "unknown"


def process_blast_file(filepath):
    """
    Read one BLAST TSV file and return summary outputs for that file.

    Parameters
    ----------
    filepath : str
        Path to one BLAST TSV file.

    Returns
    -------
    tuple
        (sample_base, top5, heatmap)

        sample_base : str | None
            The sample base name. None if the file cannot be processed.

        top5 : pandas.DataFrame | None
            A DataFrame containing the top 5 species for this sample.
            None if input could not be processed. May be empty.

        heatmap : pandas.DataFrame | None
            A DataFrame of CMV-like intervals with columns:
              sample, start, end, pident
            Empty if no CMV-like hits were found.

    Step-by-step explanation
    ------------------------
    1) Read the TSV file into a DataFrame with fixed column names.
    2) If the file is empty, return (None, None, None).
    3) Convert pident to numeric and keep rows with pident >= 85.
    4) Create a simple "species" column from stitle.
    5) Work out the sample name from the filename and store it in the DataFrame.
    6) Build a per-sample, per-species summary:
       - total_hits = number of rows (hits)
       - avg_identity = mean of pident
       Then keep the top 5.
    7) Create a CMV-only table (for heatmap-style plotting):
       - filter rows whose stitle suggests CMV
       - compute start/end from sstart/send
       - keep sample/start/end/pident
    """
    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "stitle",
            ],
        )
    except Exception as exc:
        print(f"❌ Error reading {filepath}: {exc}")
        return None, None, None

    if df.empty:
        return None, None, None

    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df = df[df["pident"] >= 85].copy()

    df["species"] = df["stitle"].apply(clean_species_name)

    sample_base, read = parse_sample_base_and_read(filepath)
    full_sample = f"{sample_base}_{read}" if read != "unknown" else sample_base
    df["sample"] = full_sample

    # Build a summary table for "top species in this sample".
    summary = (
        df.groupby(["sample", "species"])
        .agg(total_hits=("qseqid", "count"), avg_identity=("pident", "mean"))
        .reset_index()
    )

    top5 = (
        summary.sort_values(["total_hits", "avg_identity"], ascending=[False, False])
        .head(5)
        .reset_index(drop=True)
    )

    # Build a CMV-only interval table for heatmap-style plots.
    cmv_df = df[
        df["stitle"].str.contains(
            "Human herpesvirus 5|Cytomegalovirus Merlin|NC_006273.2",
            case=False,
            na=False,
        )
    ].copy()

    if not cmv_df.empty:
        cmv_df["start"] = cmv_df[["sstart", "send"]].min(axis=1)
        cmv_df["end"] = cmv_df[["sstart", "send"]].max(axis=1)
        heatmap = cmv_df[["sample", "start", "end", "pident"]].copy()
    else:
        heatmap = pd.DataFrame()

    return sample_base, top5, heatmap


def summarize_blast_hits(blast_dir, out_dir):
    """
    Summarize all BLAST TSV files in a folder and write CSV outputs.

    Parameters
    ----------
    blast_dir : str
        Directory containing BLAST TSV files named "*_blast.tsv".

    out_dir : str
        Directory to write output CSV files to. Created if missing.

    Outputs written
    ---------------
    For each sample base:
      - <sample>_blast_top5.csv
      - <sample>_blast_heatmap.csv (only if CMV hits exist)

    Also writes:
      - blast_heatmap.csv (combined heatmap for all samples)

    Step-by-step explanation
    ------------------------
    1) Make sure the output folder exists.
    2) For each TSV file:
       a) process_blast_file() reads it and returns:
          - sample base name
          - top 5 species table
          - CMV interval table
       b) Store these per sample in dictionaries.
    3) For each sample:
       a) Combine its top-5 tables (if multiple inputs exist) and write CSV.
       b) Combine its CMV interval tables (if any) and write CSV.
    4) If any CMV interval tables existed for any sample, also write one
       combined "master" heatmap CSV.
    """
    os.makedirs(out_dir, exist_ok=True)

    summary_by_sample = defaultdict(list)
    heatmap_by_sample = defaultdict(list)

    for path in sorted(glob.glob(os.path.join(blast_dir, "*_blast.tsv"))):
        sample, top5_df, heat_df = process_blast_file(path)
        if sample is None:
            continue
        if top5_df is not None and not top5_df.empty:
            summary_by_sample[sample].append(top5_df)
        if heat_df is not None and not heat_df.empty:
            heatmap_by_sample[sample].append(heat_df)

    # Write one set of output files per sample, plus a combined heatmap.
    all_heats = []

    for sample in summary_by_sample:
        combined = pd.concat(summary_by_sample[sample], ignore_index=True)
        fn = os.path.join(out_dir, f"{sample}_blast_top5.csv")
        combined.to_csv(fn, index=False)
        print(f"✅ {fn}")

        heats = heatmap_by_sample.get(sample, [])
        if heats:
            merged = pd.concat(heats, ignore_index=True)
            fn2 = os.path.join(out_dir, f"{sample}_blast_heatmap.csv")
            merged.to_csv(fn2, index=False)
            all_heats.append(merged)
            print(f"✅ {fn2}")
        else:
            print(f"⚠️ No CMV hits for {sample}")

    if all_heats:
        master = pd.concat(all_heats, ignore_index=True)
        master_fn = os.path.join(out_dir, "blast_heatmap.csv")
        master.to_csv(master_fn, index=False)
        print(f"✅ {master_fn}")


def main():
    """
    Command-line entry point.

    This function lets you run the script like:

        python3 summarize_blast_hits.py <blast_dir> [out_dir]

    Arguments
    ---------
    blast_dir:
      Folder containing "*_blast.tsv" files.

    out_dir (optional):
      Where output CSVs will be written.
      If not given, outputs are written into blast_dir.
    """
    parser = argparse.ArgumentParser(
        description="Summarize BLAST hits and create heatmap CSVs."
    )
    parser.add_argument("blast_dir", help="Directory containing *_blast.tsv files")
    parser.add_argument(
        "out_dir",
        nargs="?",
        default=None,
        help="Where to write summaries (defaults to blast_dir).",
    )
    args = parser.parse_args()

    blast_dir = args.blast_dir
    out_dir = args.out_dir or blast_dir
    summarize_blast_hits(blast_dir, out_dir)


if __name__ == "__main__":
    main()
