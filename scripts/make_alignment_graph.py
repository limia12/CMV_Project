#!/usr/bin/env python3
"""
plot_results_3_2_alignment_from_pipeline_csv.py

Creates a single dissertation-ready figure for Results section 3.2:
  - Professional violin plots (no internal median/mean/extrema lines)
  - Overlaid jittered points (one per sample)
  - Uses log10(CMV RPM + 1), derived from MAPQ-filtered CMV alignments

Designed for your pipeline outputs and batch-run workflow:
  - Reads cmv_alignment_summary.csv (preferred) or summary.csv from EACH run directory
  - Merges multiple runs (because samples were processed in separate batches)
  - Joins to truth metadata.tsv (sample -> group)
  - De-duplicates by sample ID ('base') if a sample appears in multiple runs

Inputs:
  --metadata  metadata.tsv (from your synthetic sample generator; must include sample ID + group)
  --runs      one or more run output directories containing summary.csv / cmv_alignment_summary.csv
  --outdir    output directory

Outputs:
  Fig3_2_CMV_log10RPM_by_group.png
  Table3_2_alignment_summary_by_group.tsv
  merged_alignment_metrics_for_3_2.tsv
"""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------
# Helpers
# --------------------------
def pick_col(df: pd.DataFrame, candidates) -> str:
    """Return the first matching column name from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of {candidates} found.\n"
        f"Available columns (first 80): {list(df.columns)[:80]}"
    )


def read_one_run(run_dir: Path) -> pd.DataFrame:
    """
    Read one run directory; prefer cmv_alignment_summary.csv else summary.csv.
    These are produced by your pipeline GUI write_cmv_alignment_metrics_csv().
    """
    candidates = [
        run_dir / "cmv_alignment_summary.csv",
        run_dir / "summary.csv",
    ]
    for p in candidates:
        if p.is_file() and p.stat().st_size > 0:
            df = pd.read_csv(p)
            df["_source_run"] = str(run_dir)
            df["_source_file"] = str(p)
            return df
    raise FileNotFoundError(
        f"No cmv_alignment_summary.csv or summary.csv found in: {run_dir}"
    )


def build_group_order(groups_present):
    """Preferred ordering based on your synthetic truth labels."""
    preferred = [
        "NEG_HUMAN_ONLY",
        "NEG_OTHER_VIRUS",
        "NEG_CONTAM",
        "POS_ROC_LADDER",
        "POS_CMV_PLUS_OTHER",
        "POS_MIXED_STRAINS",
        "POS_MUT_SPIKE",
    ]
    present = [g for g in preferred if g in groups_present]
    extras = [g for g in sorted(groups_present) if g not in present]
    return present + extras


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Create Results 3.2 violin plot of log10(CMV RPM + 1) by truth group."
    )
    ap.add_argument("--metadata", required=True, help="metadata.tsv containing sample and group columns")
    ap.add_argument(
        "--runs",
        required=True,
        nargs="+",
        help="One or more run directories containing summary.csv / cmv_alignment_summary.csv",
    )
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--dpi", type=int, default=450)
    ap.add_argument(
        "--title",
        default="CMV signal after host read removal and CMV alignment (MAPQ≥30)",
        help="Figure title",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------
    # Load truth metadata
    # --------------------------
    meta = pd.read_csv(args.metadata, sep="\t")
    meta_sample = pick_col(meta, ["base", "sample", "sample_id", "id"])
    meta_group = pick_col(meta, ["group", "truth_group", "label_group"])
    meta = meta.rename(columns={meta_sample: "base", meta_group: "group"})[["base", "group"]].copy()

    # --------------------------
    # Load all run summaries
    # --------------------------
    frames = []
    for r in args.runs:
        rdir = Path(r)
        if not rdir.is_dir():
            print(f"⚠️ Skipping (not a directory): {rdir}", file=sys.stderr)
            continue
        try:
            frames.append(read_one_run(rdir))
        except Exception as e:
            print(f"⚠️ Skipping run dir {rdir}: {e}", file=sys.stderr)

    if not frames:
        sys.exit("❌ No usable run summary CSVs were loaded. Check your --runs paths.")

    df = pd.concat(frames, ignore_index=True)

    # --------------------------
    # Extract required columns from pipeline summary CSV
    # --------------------------
    base_col = pick_col(df, ["base", "sample", "sample_id"])
    rpm_col = pick_col(df, ["cmv_rpm_mapq", "cmv_rpm_mapq30", "cmv_rpm"])

    df = df.rename(columns={base_col: "base", rpm_col: "cmv_rpm_mapq"}).copy()
    df["cmv_rpm_mapq"] = pd.to_numeric(df["cmv_rpm_mapq"], errors="coerce").fillna(0)

    # De-duplicate by base in case samples appear across multiple runs
    dup = int(df["base"].duplicated().sum())
    if dup > 0:
        print(
            f"⚠️ Found {dup} duplicate sample IDs across runs; keeping first occurrence per sample.",
            file=sys.stderr,
        )
        df = df.drop_duplicates(subset=["base"], keep="first")

    # --------------------------
    # Join metrics to truth groups
    # --------------------------
    merged = df.merge(meta, on="base", how="inner")
    if merged.empty:
        sys.exit(
            "❌ Join produced 0 rows. Check that metadata 'base' IDs match summary.csv 'base' IDs."
        )

    merged["log10_cmv_rpm_plus1"] = np.log10(merged["cmv_rpm_mapq"] + 1)

    merged_out = outdir / "merged_alignment_metrics_for_3_2.tsv"
    merged.to_csv(merged_out, sep="\t", index=False)

    # --------------------------
    # Build group order and summary table
    # --------------------------
    order = build_group_order(set(merged["group"].astype(str).unique()))
    merged["group"] = pd.Categorical(merged["group"], categories=order, ordered=True)
    merged = merged.sort_values("group")

    rows = []
    for g in order:
        sub = merged.loc[merged["group"] == g]
        if sub.empty:
            continue
        y = sub["log10_cmv_rpm_plus1"].astype(float).values
        pct_zero = float((sub["cmv_rpm_mapq"].values <= 0).mean() * 100.0)
        q25, q75 = np.nanpercentile(y, 25), np.nanpercentile(y, 75)
        rows.append(
            {
                "group": g,
                "n": int(len(sub)),
                "median_log10_rpm_plus1": float(np.nanmedian(y)),
                "iqr_log10_rpm_plus1": float(q75 - q25),
                "min_log10_rpm_plus1": float(np.nanmin(y)),
                "max_log10_rpm_plus1": float(np.nanmax(y)),
                "percent_zero_signal": pct_zero,
            }
        )

    summary = pd.DataFrame(rows)
    summary_path = outdir / "Table3_2_alignment_summary_by_group.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    # --------------------------
    # Plot: professional violins + jittered points (no internal lines)
    # --------------------------
    fig = plt.figure(figsize=(12, 6.5))
    ax = plt.gca()

    labels = [g for g in order if g in set(merged["group"].astype(str).unique())]
    data = [
        merged.loc[merged["group"] == g, "log10_cmv_rpm_plus1"].astype(float).values
        for g in labels
    ]

    vp = ax.violinplot(
        data,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        widths=0.85,
    )

    # Style violins: light fill, clear outline
    for body in vp["bodies"]:
        body.set_facecolor("#D0D0D0")   # light grey
        body.set_edgecolor("black")
        body.set_linewidth(0.8)
        body.set_alpha(0.7)

    # Overlay jittered points
    rng = np.random.default_rng(1)
    for i, g in enumerate(labels, start=1):
        y = merged.loc[merged["group"] == g, "log10_cmv_rpm_plus1"].astype(float).values
        x = rng.normal(loc=i, scale=0.06, size=len(y))
        ax.scatter(
            x,
            y,
            s=18,
            alpha=0.75,
            color="black",
            linewidths=0,
        )

    # Axes styling
    ax.set_title(args.title)
    ax.set_xlabel("Truth group")
    ax.set_ylabel("log10(CMV RPM + 1) (MAPQ-filtered CMV alignments)")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=35, ha="right")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    fig_path = outdir / "Fig3_2_CMV_log10RPM_by_group.png"
    plt.savefig(fig_path, dpi=args.dpi)
    plt.close(fig)

    print("✅ Done")
    print(f"Saved figure:  {fig_path}")
    print(f"Saved table:   {summary_path}")
    print(f"Saved merged:  {merged_out}")


if __name__ == "__main__":
    main()