#!/usr/bin/env python3

import os
import glob
import pandas as pd
from collections import defaultdict

def clean_species_name(stitle):
    return stitle.split(',')[0].strip()

def parse_sample_base_and_read(filename):
    basename = os.path.basename(filename)
    sample_part = basename.split("_blast")[0]
    if "_R1" in sample_part:
        return sample_part.replace("_R1", ""), "R1"
    elif "_R2" in sample_part:
        return sample_part.replace("_R2", ""), "R2"
    else:
        return sample_part, "unknown"

def process_blast_file(filepath):
    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=[
                "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
            ]
        )
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        return None, None, None

    if df.empty:
        return None, None, None

    df = df[df["pident"] >= 85].copy()
    df["species"] = df["stitle"].apply(clean_species_name)

    sample_base, read = parse_sample_base_and_read(filepath)
    full_sample = f"{sample_base}_{read}"
    df["sample"] = full_sample

    # === Top 5 Species Summary ===
    summary = (
        df.groupby(["sample", "species"])
          .agg(total_hits=("qseqid", "count"), avg_identity=("pident", "mean"))
          .reset_index()
    )
    top5 = (
        summary.sort_values(by=["total_hits", "avg_identity"], ascending=[False, False])
               .head(5)
               .reset_index(drop=True)
    )

    # === CMV-Only Heatmap ===
    cmv_df = df[df["stitle"].str.contains("Human herpesvirus 5|Cytomegalovirus Merlin|NC_006273.2", case=False, na=False)]
    if not cmv_df.empty:
        cmv_df["start"] = cmv_df[["sstart", "send"]].min(axis=1)
        cmv_df["end"] = cmv_df[["sstart", "send"]].max(axis=1)
        heatmap = cmv_df[["sample", "start", "end", "pident"]].copy()
    else:
        heatmap = pd.DataFrame()

    return sample_base, top5, heatmap

def summarize_blast_hits(blast_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # Store by sample base name
    summary_by_sample = defaultdict(list)
    heatmap_by_sample = defaultdict(list)

    for filepath in glob.glob(os.path.join(blast_dir, "*_blast.tsv")):
        sample_base, top5_df, heatmap_df = process_blast_file(filepath)
        if sample_base is None:
            continue
        summary_by_sample[sample_base].append(top5_df)
        if heatmap_df is not None and not heatmap_df.empty:
            heatmap_by_sample[sample_base].append(heatmap_df)

    # Write one summary and one heatmap file per sample_base
    for sample_base in summary_by_sample:
        combined_top5 = pd.concat(summary_by_sample[sample_base], ignore_index=True)
        combined_top5_path = os.path.join(out_dir, f"{sample_base}_blast_top5.csv")
        combined_top5.to_csv(combined_top5_path, index=False)
        print(f"✅ Saved top 5 summary: {combined_top5_path}")

        if sample_base in heatmap_by_sample:
            combined_heatmap = pd.concat(heatmap_by_sample[sample_base], ignore_index=True)
            heatmap_path = os.path.join(out_dir, f"{sample_base}_blast_heatmap.csv")
            combined_heatmap.to_csv(heatmap_path, index=False)
            print(f"✅ Saved heatmap data: {heatmap_path}")
        else:
            print(f"⚠️ No CMV hits for {sample_base}, no heatmap generated.")

