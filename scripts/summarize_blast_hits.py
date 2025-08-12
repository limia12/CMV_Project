#!/usr/bin/env python3
import os
import sys
import glob
import pandas as pd
from collections import defaultdict
import re
import argparse

def clean_species_name(stitle):
    return stitle.split(',')[0].strip()

def parse_sample_base_and_read(filename):
    basename = os.path.basename(filename)
    sample_part = basename.split("_blast")[0]
    if re.search(r"_R[12]$", sample_part):
        return sample_part.rsplit("_", 1)[0], sample_part.rsplit("_", 1)[1]
    elif "_R1" in sample_part or "_R2" in sample_part:
        return sample_part.replace("_R1", "").replace("_R2", ""), "unknown"
    else:
        return sample_part, "unknown"

def process_blast_file(filepath):
    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=[
                "qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore","stitle"
            ]
        )
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        return None, None, None

    if df.empty:
        return None, None, None

    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df = df[df["pident"] >= 85].copy()
    df["species"] = df["stitle"].apply(clean_species_name)

    sample_base, read = parse_sample_base_and_read(filepath)
    full_sample = f"{sample_base}_{read}" if read != "unknown" else sample_base
    df["sample"] = full_sample

    # Top-5 species summary
    summary = (
        df.groupby(["sample","species"])
          .agg(total_hits=("qseqid","count"), avg_identity=("pident","mean"))
          .reset_index()
    )
    top5 = (
        summary.sort_values(["total_hits","avg_identity"], ascending=[False,False])
               .head(5)
               .reset_index(drop=True)
    )

    # CMV‐specific intervals for heatmap
    cmv_df = df[df["stitle"].str.contains(
        "Human herpesvirus 5|Cytomegalovirus Merlin|NC_006273.2",
        case=False, na=False
    )].copy()
    if not cmv_df.empty:
        cmv_df["start"] = cmv_df[["sstart","send"]].min(axis=1)
        cmv_df["end"]   = cmv_df[["sstart","send"]].max(axis=1)
        heatmap = cmv_df[["sample","start","end","pident"]].copy()
    else:
        heatmap = pd.DataFrame()

    return sample_base, top5, heatmap

def summarize_blast_hits(blast_dir, out_dir):
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

    # Write per-sample summary & heatmap, and aggregate a unified file
    all_heats = []
    for sample in summary_by_sample:
        # Top-5 CSV
        combined = pd.concat(summary_by_sample[sample], ignore_index=True)
        fn = os.path.join(out_dir, f"{sample}_blast_top5.csv")
        combined.to_csv(fn, index=False)
        print(f"✅ {fn}")

        # Interval heatmap CSV
        heats = heatmap_by_sample.get(sample, [])
        if heats:
            merged = pd.concat(heats, ignore_index=True)
            fn2 = os.path.join(out_dir, f"{sample}_blast_heatmap.csv")
            merged.to_csv(fn2, index=False)
            all_heats.append(merged)
            print(f"✅ {fn2}")
        else:
            print(f"⚠️ No CMV hits for {sample}")

    # Unified master heatmap
    if all_heats:
        master = pd.concat(all_heats, ignore_index=True)
        master_fn = os.path.join(out_dir, "blast_heatmap.csv")
        master.to_csv(master_fn, index=False)
        print(f"✅ {master_fn}")

def main():
    p = argparse.ArgumentParser(description="Summarize BLAST hits & build heatmap CSVs")
    p.add_argument("blast_dir", help="Directory containing *_blast.tsv files")
    p.add_argument("out_dir", nargs="?", default=None,
                   help="Where to write summaries (defaults to blast_dir)")
    args = p.parse_args()

    blast_dir = args.blast_dir
    out_dir   = args.out_dir or blast_dir
    summarize_blast_hits(blast_dir, out_dir)

if __name__ == "__main__":
    main()
