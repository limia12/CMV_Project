#!/usr/bin/env python3
"""
make_coverage_graph.py

Creates a single dissertation-ready figure for Results section 3.3:
  - Representative conserved-gene coverage (IE1, UL54, UL97, UL83)
  - One true positive vs one contamination vs one negative
  - Designed for batch-run workflow (samples processed across multiple run folders)

Inputs:
  --metadata  metadata.tsv (truth table; must include sample ID + group)
  --runs      one or more run output directories (searched recursively for *_gene.csv)
  --outdir    output directory

Outputs:
  Fig3_3_representative_gene_coverage.png
  Table3_3_selected_representative_samples.tsv
  merged_gene_coverage_for_3_3.tsv

Notes:
  - Uses the SAME canonical sample ID stripping logic as your Streamlit app (app.py),
    so joins should work even if files contain suffixes like _R1, _sorted, .bam, etc.
"""

import argparse
from pathlib import Path
import sys
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------
# Canonical sample naming (copied from app.py)
# --------------------------
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
    r"(?:\.bam)$",
]


def canonical_sample(name: str) -> str:
    """Strip common pipeline/file suffixes repeatedly until stable."""
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


# --------------------------
# Helpers
# --------------------------
def pick_col(df: pd.DataFrame, candidates) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of {candidates} found.\nAvailable columns (first 80): {list(df.columns)[:80]}"
    )


def find_gene_files(run_dir: Path):
    """Recursively find *_gene.csv within a run directory."""
    return sorted(run_dir.rglob("*_gene.csv"))


def read_gene_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["_source_file"] = str(path)
    df["_source_run"] = str(path.parent)
    return df


def choose_group_names(groups_present):
    """Pick POS, CONTAM, NEG groups with sensible defaults."""
    groups_present = [str(g) for g in groups_present]

    # contamination
    contam = None
    if "NEG_CONTAM" in groups_present:
        contam = "NEG_CONTAM"
    else:
        for g in groups_present:
            if "CONTAM" in g.upper():
                contam = g
                break

    # negative
    neg = "NEG_HUMAN_ONLY" if "NEG_HUMAN_ONLY" in groups_present else None
    if neg is None:
        negs = [g for g in groups_present if g.startswith("NEG_")]
        neg = negs[0] if negs else None

    # positive
    poss = [g for g in groups_present if g.startswith("POS_")]
    pos = poss[0] if poss else None

    return pos, contam, neg


def pick_representative_sample(df_group: pd.DataFrame, depth_col: str) -> str:
    """
    Pick a representative sample ID from a group.
    Uses the median overall mean depth across the conserved genes.
    """
    # overall per-sample depth summary across the conserved genes
    per_sample = (
        df_group.groupby("base")[depth_col]
        .mean()
        .reset_index()
        .sort_values(depth_col)
    )
    if per_sample.empty:
        return None
    # median row
    return per_sample.iloc[len(per_sample) // 2]["base"]


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Results 3.3: Representative conserved-gene coverage plot from *_gene.csv files."
    )
    ap.add_argument("--metadata", required=True, help="metadata.tsv containing sample and group columns")
    ap.add_argument("--runs", required=True, nargs="+", help="One or more run output directories")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--dpi", type=int, default=450)
    ap.add_argument(
        "--genes",
        nargs="+",
        default=["IE1", "UL54", "UL97", "UL83"],
        help="Genes to plot (default: IE1 UL54 UL97 UL83)",
    )
    ap.add_argument(
        "--log10",
        action="store_true",
        help="Plot log10(mean_depth + 1) instead of raw mean depth (recommended for readability).",
    )
    ap.add_argument(
        "--title",
        default="Representative conserved-gene coverage patterns (report output)",
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
    meta = meta.rename(columns={meta_sample: "base_raw", meta_group: "group"})[["base_raw", "group"]].copy()
    meta["base"] = meta["base_raw"].astype(str).apply(canonical_sample)

    # --------------------------
    # Load all *_gene.csv across runs
    # --------------------------
    gene_frames = []
    searched = 0
    for r in args.runs:
        rdir = Path(r)
        if not rdir.is_dir():
            print(f"⚠️ Skipping (not a directory): {rdir}", file=sys.stderr)
            continue
        files = find_gene_files(rdir)
        searched += len(files)
        for f in files:
            try:
                if f.stat().st_size == 0:
                    continue
                gene_frames.append(read_gene_file(f))
            except Exception as e:
                print(f"⚠️ Failed reading {f}: {e}", file=sys.stderr)

    if not gene_frames:
        sys.exit(
            "❌ No *_gene.csv files were loaded.\n"
            "Check that your --runs directories contain gene coverage outputs (pattern: *_gene.csv)."
        )

    genes = pd.concat(gene_frames, ignore_index=True)

    # Identify columns in gene table
    sample_col = pick_col(genes, ["sample", "base", "sample_id"])
    gene_col = pick_col(genes, ["gene", "region", "target"])
    depth_col = pick_col(genes, ["mean_depth", "avg_depth", "average_depth", "depth_mean", "mean"])

    genes = genes.rename(columns={sample_col: "sample", gene_col: "gene", depth_col: "mean_depth"}).copy()
    genes["base"] = genes["sample"].astype(str).apply(canonical_sample)
    genes["mean_depth"] = pd.to_numeric(genes["mean_depth"], errors="coerce")

    # Keep only the genes we want
    genes = genes[genes["gene"].astype(str).isin(args.genes)].copy()
    if genes.empty:
        sys.exit(
            f"❌ After filtering, no rows matched genes: {args.genes}\n"
            f"Check what gene names your *_gene.csv uses (e.g. gpUL55 vs UL55)."
        )

    # De-duplicate: if the same base+gene appears across runs, keep the first non-null mean_depth
    genes = genes.sort_values(["base", "gene"]).drop_duplicates(subset=["base", "gene"], keep="first")

    # --------------------------
    # Join to truth groups (CRITICAL)
    # --------------------------
    merged = genes.merge(meta[["base", "group"]], on="base", how="inner")
    if merged.empty:
        # debug: show what IDs exist on each side
        meta_ids = set(meta["base"].dropna().astype(str).unique())
        gene_ids = set(genes["base"].dropna().astype(str).unique())
        inter = sorted(meta_ids.intersection(gene_ids))

        debug_path = outdir / "DEBUG_base_id_sets.txt"
        with open(debug_path, "w") as fh:
            fh.write("=== Example metadata base IDs (first 30) ===\n")
            fh.write("\n".join(sorted(list(meta_ids))[:30]) + "\n\n")
            fh.write("=== Example gene-table base IDs (first 30) ===\n")
            fh.write("\n".join(sorted(list(gene_ids))[:30]) + "\n\n")
            fh.write(f"=== Intersection size: {len(inter)} ===\n")
            fh.write("\n".join(inter[:50]) + "\n")

        sys.exit(
            "❌ Join produced 0 rows.\n"
            "This means the sample IDs still don't match after canonicalisation.\n"
            f"I wrote a debug file you can inspect:\n  {debug_path}\n"
            "Common causes:\n"
            "  - metadata sample IDs include extra prefixes not present in filenames (or vice versa)\n"
            "  - gene files use a different identifier than metadata (e.g. full path, different naming convention)\n"
        )

    merged.to_csv(outdir / "merged_gene_coverage_for_3_3.tsv", sep="\t", index=False)

    # --------------------------
    # Choose representative groups + samples
    # --------------------------
    groups_present = sorted(merged["group"].dropna().astype(str).unique())
    pos_group, contam_group, neg_group = choose_group_names(groups_present)

    if pos_group is None or contam_group is None or neg_group is None:
        sys.exit(
            "❌ Could not identify POS / CONTAM / NEG groups automatically.\n"
            f"Groups present were:\n  {groups_present}\n"
            "If your labels differ, rename them in metadata.tsv (recommended), "
            "or edit choose_group_names() in this script."
        )

    pos_sample = pick_representative_sample(merged[merged["group"] == pos_group], "mean_depth")
    contam_sample = pick_representative_sample(merged[merged["group"] == contam_group], "mean_depth")
    neg_sample = pick_representative_sample(merged[merged["group"] == neg_group], "mean_depth")

    if not all([pos_sample, contam_sample, neg_sample]):
        sys.exit("❌ Could not pick one representative sample per group (insufficient data).")

    chosen = pd.DataFrame(
        [
            {"role": "True positive (representative)", "group": pos_group, "base": pos_sample},
            {"role": "Contamination (representative)", "group": contam_group, "base": contam_sample},
            {"role": "Negative (representative)", "group": neg_group, "base": neg_sample},
        ]
    )
    chosen.to_csv(outdir / "Table3_3_selected_representative_samples.tsv", sep="\t", index=False)

    # Subset merged to just the chosen samples
    keep_bases = [pos_sample, contam_sample, neg_sample]
    sub = merged[merged["base"].isin(keep_bases)].copy()

    # Make plotting values
    if args.log10:
        sub["plot_depth"] = np.log10(sub["mean_depth"].fillna(0) + 1)
        ylab = "log10(mean depth + 1)"
    else:
        sub["plot_depth"] = sub["mean_depth"]
        ylab = "Mean depth"

    # Ensure gene order as requested
    gene_order = args.genes
    sub["gene"] = pd.Categorical(sub["gene"], categories=gene_order, ordered=True)
    sub = sub.sort_values("gene")

    # Pivot to gene x sample matrix
    mat = (
        sub.pivot_table(index="gene", columns="base", values="plot_depth", aggfunc="mean")
        .reindex(gene_order)
    )

    # --------------------------
    # Plot: grouped bar chart (professional, simple)
    # --------------------------
    fig = plt.figure(figsize=(10.8, 5.8))
    ax = plt.gca()

    x = np.arange(len(gene_order))
    bases = keep_bases
    width = 0.26

    for i, b in enumerate(bases):
        y = mat[b].values if b in mat.columns else np.zeros(len(gene_order))
        ax.bar(x + (i - 1) * width, y, width=width, label=b)

    ax.set_title(args.title)
    ax.set_xlabel("Conserved CMV genes")
    ax.set_ylabel(ylab)
    ax.set_xticks(x)
    ax.set_xticklabels(gene_order)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Put legend outside for clarity
    ax.legend(title="Representative sample (base ID)", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()

    fig_path = outdir / "Fig3_3_representative_gene_coverage.png"
    plt.savefig(fig_path, dpi=args.dpi)
    plt.close(fig)

    print("✅ Done")
    print(f"Selected groups: POS={pos_group} | CONTAM={contam_group} | NEG={neg_group}")
    print(f"Selected samples: POS={pos_sample} | CONTAM={contam_sample} | NEG={neg_sample}")
    print(f"Saved figure: {fig_path}")
    print(f"Saved table:  {outdir / 'Table3_3_selected_representative_samples.tsv'}")
    print(f"Saved merged: {outdir / 'merged_gene_coverage_for_3_3.tsv'}")


if __name__ == "__main__":
    main()