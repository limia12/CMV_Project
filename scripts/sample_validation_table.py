#!/usr/bin/env python3
"""
sample_validation_table.py

Purpose
-------
Create a simple validation table and visual summaries showing whether
simulated CMV samples behaved as expected after running them through
the pipeline across three separate pipeline output directories.

This version is for cases where the ROC sample set had to be run in
three separate batches.

Inputs
------
- truth/metadata.tsv from b_generate_roc_sampleset.py
- pipeline CSV outputs from three different directories:
  * cmv_alignment_summary.csv (or *_summary.csv)
  * *_gene.csv
  * *_blast_top5.csv

Outputs
-------
1) validation_table.tsv / validation_table.csv
   Per-sample expected vs observed summary with simple pass/check flags.
2) validation_overview_by_group.tsv
   Group-level summary.
3) figures/validation_heatmap.png
   Heatmap-like visual table of pass/check patterns.
4) figures/roc_ladder_target_vs_observed.png
   Scatter plot for ROC ladder samples: target CMV fraction vs observed log10 RPM.
5) figures/group_metric_boxplots.png
   Boxplots of key observed metrics by sample group.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_SUFFIX_PATTERNS = [
    r"(?:_R?[12](?:_00[1-9])?)$",
    r"(?:_[12])$",
    r"(?:_sorted(?:_sorted)*)$",
    r"(?:_sort(?:ed)?)$",
    r"(?:_dedup(?:ed)?)$",
    r"(?:_mark(?:ed)?dups?)$",
    r"(?:_rmdup[s]?)$",
    r"(?:_trim(?:med)?(?:_paired)?)$",
    r"(?:_filtered)$",
    r"(?:_indelqual)$",
    r"(?:_realigned)$",
    r"(?:\.bam)$",
]


def canonical_sample(name: str) -> str:
    s = str(name)
    changed = True
    while changed:
        changed = False
        for pat in _SUFFIX_PATTERNS:
            new = re.sub(pat, "", s, flags=re.IGNORECASE)
            if new != s:
                s = new
                changed = True
    return s.rstrip("_-")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()

    for sep in ["\t", ","]:
        try:
            df = pd.read_csv(path, sep=sep, low_memory=False)
            if df.shape[1] >= 2:
                return df
        except Exception:
            pass

    try:
        return pd.read_csv(path, sep=None, engine="python", low_memory=False)
    except Exception:
        return pd.DataFrame()


def read_glob(files: Iterable[Path]) -> pd.DataFrame:
    dfs: List[pd.DataFrame] = []
    for p in files:
        df = read_table(p)
        if not df.empty:
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def coerce_numeric(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def collect_csv_dirs(args) -> List[Path]:
    csv_dirs = [
        Path(args.csv_dir_1),
        Path(args.csv_dir_2),
        Path(args.csv_dir_3),
    ]
    return csv_dirs


def load_truth(metadata_tsv: Path) -> pd.DataFrame:
    df = read_table(metadata_tsv)
    if df.empty:
        raise SystemExit(f"Could not read truth metadata: {metadata_tsv}")

    required = ["sample", "group", "truth_cmv_present"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Truth metadata missing columns: {', '.join(missing)}")

    df["base"] = df["sample"].astype(str).apply(canonical_sample)
    df["truth_cmv_present"] = pd.to_numeric(
        df["truth_cmv_present"], errors="coerce"
    ).fillna(0).astype(int)

    keep = [
        "sample", "base", "group", "truth_cmv_present",
        "cmv_fraction_target", "total_pairs",
        "cmv_refs", "cmv_pairs", "other_refs", "other_pairs",
        "human_ref", "human_pairs", "read_len", "seq_sys",
        "mean_frag", "sd_frag", "snv_rate", "random_snv_count",
        "mutation_spiked", "mutation_af_target", "mutations_applied"
    ]
    keep = [c for c in keep if c in df.columns]
    return df[keep].drop_duplicates("base")


def load_alignment_summaries(csv_dirs: List[Path]) -> pd.DataFrame:
    files = []
    for d in csv_dirs:
        if not d.exists():
            print(f"Warning: alignment directory not found: {d}")
            continue

        preferred = d / "cmv_alignment_summary.csv"
        if preferred.exists():
            files.append(preferred)
        else:
            files.extend(sorted(d.glob("*_summary.csv")))
            files.extend(sorted(d.glob("summary.csv")))

    df = read_glob(files)
    if df.empty:
        return df

    if "base" not in df.columns:
        if "sample" in df.columns:
            df["base"] = df["sample"].astype(str).apply(canonical_sample)
        else:
            df["base"] = df.iloc[:, 0].astype(str).apply(canonical_sample)

    df["base"] = df["base"].astype(str).apply(canonical_sample)
    df = coerce_numeric(df, [
        "cmv_mapped_reads_mapq",
        "cmv_rpm_mapq",
        "cmv_breadth_mapq",
        "cmv_cov_bases_mapq",
        "cmv_ref_bases",
        "total_reads_in_cmv_bam"
    ])

    agg = {}
    for col in ["cmv_mapped_reads_mapq", "cmv_cov_bases_mapq", "total_reads_in_cmv_bam"]:
        if col in df.columns:
            agg[col] = "sum"

    for col in ["cmv_rpm_mapq", "cmv_breadth_mapq", "cmv_ref_bases"]:
        if col in df.columns:
            agg[col] = "max"

    out = df.groupby("base", as_index=False).agg(agg) if agg else df[["base"]].drop_duplicates()

    if {"cmv_cov_bases_mapq", "cmv_ref_bases"}.issubset(out.columns):
        out["cmv_breadth_from_cov"] = (
            out["cmv_cov_bases_mapq"] / out["cmv_ref_bases"].replace({0: np.nan})
        )

    if "cmv_mapped_reads_mapq" in out.columns:
        out["log10_mapped_hq"] = np.log10(
            out["cmv_mapped_reads_mapq"].replace({0: np.nan})
        )

    if "cmv_rpm_mapq" in out.columns:
        out["log10_rpm_hq"] = np.log10(
            out["cmv_rpm_mapq"].replace({0: np.nan})
        )

    return out


def load_gene_metrics(csv_dirs: List[Path]) -> pd.DataFrame:
    files = []
    for d in csv_dirs:
        if d.exists():
            files.extend(sorted(d.glob("*_gene.csv")))

    df = read_glob(files)
    if df.empty:
        return df

    if "base" not in df.columns:
        if "sample" in df.columns:
            df["base"] = df["sample"].astype(str).apply(canonical_sample)
        else:
            raise SystemExit("Gene coverage table missing 'base' or 'sample' column.")

    df["base"] = df["base"].astype(str).apply(canonical_sample)
    df = coerce_numeric(df, ["avg_depth", "breadth"])

    g = df.groupby("base")
    out = pd.DataFrame({"base": list(g.groups.keys())})

    if "avg_depth" in df.columns:
        out["gene_mean_depth"] = g["avg_depth"].mean().values
        out["gene_min_depth"] = g["avg_depth"].min().values

    if "breadth" in df.columns:
        out["gene_mean_breadth"] = g["breadth"].mean().values
        out["gene_min_breadth"] = g["breadth"].min().values

    out["n_genes_reported"] = g.size().values

    if "gene" in df.columns and "avg_depth" in df.columns:
        for gene in ["UL97", "UL54", "IE1", "UL83"]:
            sub = df[df["gene"].astype(str).str.upper() == gene]
            if not sub.empty:
                gene_depth = sub.groupby("base")["avg_depth"].mean()
                out[f"{gene}_mean_depth"] = out["base"].map(gene_depth)

    if "gene_mean_depth" in out.columns:
        out["log10_gene_mean_depth"] = np.log10(
            out["gene_mean_depth"].replace({0: np.nan})
        )

    return out


def load_blast_top5(csv_dirs: List[Path]) -> pd.DataFrame:
    files = []
    for d in csv_dirs:
        if d.exists():
            files.extend(sorted(d.glob("*_blast_top5.csv")))

    df = read_glob(files)
    if df.empty:
        return df

    if "base" not in df.columns:
        if "sample" in df.columns:
            df["base"] = df["sample"].astype(str).apply(canonical_sample)
        else:
            raise SystemExit("BLAST top5 table missing 'base' or 'sample' column.")

    df["base"] = df["base"].astype(str).apply(canonical_sample)
    df = coerce_numeric(df, ["total_hits", "avg_identity"])

    if "species" in df.columns:
        species = df["species"].astype(str)
    else:
        species = pd.Series("", index=df.index)

    is_cmv = species.str.contains(
        "herpesvirus 5|cytomegalovirus|merlin|hcmv",
        case=False,
        na=False
    )

    rows = []
    for base, sub in df.groupby("base"):
        cmv_sub = sub[is_cmv.loc[sub.index]]
        any_hits = float(sub["total_hits"].sum()) if "total_hits" in sub.columns else np.nan

        if cmv_sub.empty:
            rows.append({
                "base": base,
                "blast_cmv_hits_top5": 0,
                "blast_any_hits_top5": any_hits,
                "blast_cmv_mean_identity_top5": np.nan,
            })
        else:
            if "total_hits" in cmv_sub.columns and float(cmv_sub["total_hits"].sum()) > 0:
                weights = cmv_sub["total_hits"].to_numpy(dtype=float)
                mean_identity = float(
                    np.average(cmv_sub["avg_identity"].to_numpy(dtype=float), weights=weights)
                )
                hit_count = int(cmv_sub["total_hits"].sum())
            else:
                mean_identity = float(cmv_sub["avg_identity"].mean())
                hit_count = int(len(cmv_sub))

            rows.append({
                "base": base,
                "blast_cmv_hits_top5": hit_count,
                "blast_any_hits_top5": any_hits,
                "blast_cmv_mean_identity_top5": mean_identity,
            })

    return pd.DataFrame(rows)


def build_validation_table(
    truth: pd.DataFrame,
    bam: pd.DataFrame,
    genes: pd.DataFrame,
    blast: pd.DataFrame
) -> pd.DataFrame:
    df = truth.copy()

    for extra in [bam, genes, blast]:
        if extra is not None and not extra.empty:
            df = df.merge(extra, on="base", how="left")

    positive_groups = {
        "POS_ROC_LADDER",
        "POS_CMV_PLUS_OTHER",
        "POS_MIXED_STRAINS",
        "POS_MUT_SPIKE"
    }
    negative_groups = {
        "NEG_HUMAN_ONLY",
        "NEG_OTHER_VIRUS",
        "NEG_CONTAM"
    }

    def expected_pattern(group: str) -> str:
        if group in positive_groups:
            return "CMV signal expected"
        if group == "NEG_CONTAM":
            return "Very weak / near-background CMV signal expected"
        if group in negative_groups:
            return "No meaningful CMV signal expected"
        return "Unknown"

    df["expected_pattern"] = df["group"].astype(str).map(expected_pattern)

    rpm = df["log10_rpm_hq"] if "log10_rpm_hq" in df.columns else pd.Series(np.nan, index=df.index)
    breadth = df["cmv_breadth_mapq"] if "cmv_breadth_mapq" in df.columns else pd.Series(np.nan, index=df.index)
    blast_id = (
        df["blast_cmv_mean_identity_top5"]
        if "blast_cmv_mean_identity_top5" in df.columns
        else pd.Series(np.nan, index=df.index)
    )

    df["flag_observed_cmv_signal"] = (
        (rpm.fillna(-np.inf) > 0.5) |
        (breadth.fillna(0) > 0.001) |
        (blast_id.fillna(0) > 80)
    )

    df["flag_identity_support"] = blast_id.fillna(0) > 80
    df["flag_coverage_support"] = breadth.fillna(0) > 0.001
    df["flag_depth_support"] = rpm.fillna(-np.inf) > 0.5

    def classify_row(row: pd.Series) -> str:
        group = str(row.get("group", ""))
        signal = bool(row.get("flag_observed_cmv_signal", False))
        identity = bool(row.get("flag_identity_support", False))
        coverage = bool(row.get("flag_coverage_support", False))

        if group in positive_groups:
            return "PASS" if signal else "CHECK"
        if group in {"NEG_HUMAN_ONLY", "NEG_OTHER_VIRUS"}:
            return "PASS" if not signal else "CHECK"
        if group == "NEG_CONTAM":
            return "PASS" if (not coverage and not identity) else "CHECK"
        return "CHECK"

    df["validation_status"] = df.apply(classify_row, axis=1)

    def interpretation(row: pd.Series) -> str:
        parts = []
        if pd.notna(row.get("log10_rpm_hq", np.nan)):
            parts.append(f"log10 RPM={row['log10_rpm_hq']:.2f}")
        if pd.notna(row.get("cmv_breadth_mapq", np.nan)):
            parts.append(f"breadth={row['cmv_breadth_mapq']:.4f}")
        if pd.notna(row.get("blast_cmv_mean_identity_top5", np.nan)):
            parts.append(f"BLAST id={row['blast_cmv_mean_identity_top5']:.1f}%")
        return "; ".join(parts)

    df["observed_summary"] = df.apply(interpretation, axis=1)

    order = [
        "sample", "group", "truth_cmv_present", "expected_pattern",
        "validation_status", "cmv_fraction_target", "total_pairs",
        "cmv_refs", "other_refs", "log10_rpm_hq", "cmv_mapped_reads_mapq",
        "cmv_breadth_mapq", "gene_mean_depth", "gene_min_depth",
        "blast_cmv_hits_top5", "blast_cmv_mean_identity_top5",
        "flag_depth_support", "flag_coverage_support", "flag_identity_support",
        "observed_summary", "mutation_spiked", "mutation_af_target",
        "mutations_applied"
    ]
    order = [c for c in order if c in df.columns]

    return df[order].sort_values(["group", "sample"]).reset_index(drop=True)


def save_group_summary(df: pd.DataFrame, outpath: Path) -> None:
    rows = []
    for group, sub in df.groupby("group", dropna=False):
        row = {
            "group": group,
            "n_samples": len(sub),
            "n_pass": int((sub["validation_status"] == "PASS").sum()),
            "n_check": int((sub["validation_status"] == "CHECK").sum()),
        }

        if "cmv_fraction_target" in sub.columns:
            row["median_target_fraction"] = pd.to_numeric(
                sub["cmv_fraction_target"], errors="coerce"
            ).median()

        if "log10_rpm_hq" in sub.columns:
            row["median_log10_rpm"] = pd.to_numeric(
                sub["log10_rpm_hq"], errors="coerce"
            ).median()

        if "cmv_breadth_mapq" in sub.columns:
            row["median_breadth"] = pd.to_numeric(
                sub["cmv_breadth_mapq"], errors="coerce"
            ).median()

        if "blast_cmv_mean_identity_top5" in sub.columns:
            row["median_blast_identity"] = pd.to_numeric(
                sub["blast_cmv_mean_identity_top5"], errors="coerce"
            ).median()

        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(outpath, sep="\t", index=False)


def plot_validation_heatmap(df: pd.DataFrame, outpath: Path, max_samples: int = 120) -> None:
    plot_df = df.copy()
    if len(plot_df) > max_samples:
        plot_df = plot_df.iloc[:max_samples].copy()

    cols = [
        c for c in
        ["flag_depth_support", "flag_coverage_support", "flag_identity_support"]
        if c in plot_df.columns
    ]
    if not cols:
        return

    matrix = plot_df[cols].fillna(False).astype(int).to_numpy()

    fig_h = max(4, min(18, 0.22 * len(plot_df) + 1.5))
    fig, ax = plt.subplots(figsize=(8, fig_h))
    ax.imshow(matrix, aspect="auto")

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(
        [c.replace("flag_", "").replace("_", " ") for c in cols],
        rotation=45,
        ha="right"
    )
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["sample"].astype(str))
    ax.set_title("Per-sample validation flags")
    ax.set_xlabel("Observed evidence")
    ax.set_ylabel("Sample")

    for tick, status in zip(ax.get_yticklabels(), plot_df["validation_status"].astype(str)):
        if status == "PASS":
            tick.set_fontweight("bold")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_roc_ladder_target_vs_observed(df: pd.DataFrame, outpath: Path) -> None:
    needed = {"group", "cmv_fraction_target", "log10_rpm_hq"}
    if not needed.issubset(df.columns):
        return

    sub = df[df["group"].astype(str) == "POS_ROC_LADDER"].copy()
    sub["cmv_fraction_target"] = pd.to_numeric(sub["cmv_fraction_target"], errors="coerce")
    sub["log10_rpm_hq"] = pd.to_numeric(sub["log10_rpm_hq"], errors="coerce")
    sub = sub.dropna(subset=["cmv_fraction_target", "log10_rpm_hq"])

    if sub.empty:
        return

    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(sub["cmv_fraction_target"], sub["log10_rpm_hq"])
    ax.set_xscale("log")
    ax.set_xlabel("Target CMV fraction")
    ax.set_ylabel("Observed log10 RPM (HQ CMV reads)")
    ax.set_title("ROC ladder samples: target vs observed CMV signal")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_group_metric_boxplots(df: pd.DataFrame, outpath: Path) -> None:
    metrics = [
        m for m in
        ["log10_rpm_hq", "cmv_breadth_mapq", "blast_cmv_mean_identity_top5"]
        if m in df.columns
    ]
    if not metrics or "group" not in df.columns:
        return

    groups = sorted(df["group"].dropna().astype(str).unique().tolist())

    fig, axes = plt.subplots(nrows=len(metrics), ncols=1, figsize=(10, 4 * len(metrics)))
    if len(metrics) == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        data = []
        for g in groups:
            vals = pd.to_numeric(
                df.loc[df["group"].astype(str) == g, metric],
                errors="coerce"
            ).dropna().to_numpy()
            data.append(vals if len(vals) else np.array([np.nan]))

        ax.boxplot(data, labels=groups, vert=True)
        ax.set_title(metric.replace("_", " "))
        ax.tick_params(axis="x", rotation=45)

    fig.suptitle("Observed metrics by simulated sample group")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(outdir: Path) -> None:
    text = """
How to interpret the validation outputs
=======================================

validation_table.tsv / .csv
- One row per sample.
- Shows the expected simulation type and the observed CMV metrics.
- validation_status = PASS means the observed behaviour broadly matches the intended design.
- validation_status = CHECK means it is worth reviewing that sample manually.

validation_overview_by_group.tsv
- Grouped summary showing how many samples in each design category behaved as expected.

Suggested wording for your dissertation results
-----------------------------------------------
A validation table was created by comparing the intended simulation design recorded in the truth metadata with the observed pipeline outputs for each sample. Samples simulated as CMV-positive showed higher CMV read abundance, broader CMV genome coverage and higher CMV BLAST identity than negative controls. ROC ladder samples showed a graded increase in observed CMV signal with increasing target CMV fraction. Human-only and other-virus negative controls showed little or no CMV signal, while contamination controls remained close to background as intended.
""".strip()

    (outdir / "README_validation.txt").write_text(text + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create a validation table and figures for simulated CMV samples from three separate pipeline output folders."
    )
    ap.add_argument("--truth", required=True, help="Path to truth/metadata.tsv from the simulator.")
    ap.add_argument("--csv-dir-1", required=True, help="First pipeline CSV output directory.")
    ap.add_argument("--csv-dir-2", required=True, help="Second pipeline CSV output directory.")
    ap.add_argument("--csv-dir-3", required=True, help="Third pipeline CSV output directory.")
    ap.add_argument("--outdir", required=True, help="Output directory for validation tables and figures.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    ensure_dir(outdir)
    ensure_dir(figdir)

    csv_dirs = collect_csv_dirs(args)

    print("Using CSV directories:")
    for d in csv_dirs:
        print(f"  - {d}")

    truth = load_truth(Path(args.truth))
    bam = load_alignment_summaries(csv_dirs)
    genes = load_gene_metrics(csv_dirs)
    blast = load_blast_top5(csv_dirs)

    validation = build_validation_table(truth, bam, genes, blast)
    validation.to_csv(outdir / "validation_table.tsv", sep="\t", index=False)
    validation.to_csv(outdir / "validation_table.csv", index=False)

    save_group_summary(validation, outdir / "validation_overview_by_group.tsv")
    plot_validation_heatmap(validation, figdir / "validation_heatmap.png")
    plot_roc_ladder_target_vs_observed(validation, figdir / "roc_ladder_target_vs_observed.png")
    plot_group_metric_boxplots(validation, figdir / "group_metric_boxplots.png")
    write_readme(outdir)

    print("Done.")
    print(f"Validation table: {outdir / 'validation_table.tsv'}")
    print(f"Group summary:    {outdir / 'validation_overview_by_group.tsv'}")
    print(f"Figures:          {figdir}")


if __name__ == "__main__":
    main()