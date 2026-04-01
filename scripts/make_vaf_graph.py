#!/usr/bin/env python3
"""
make_vaf_graph.py

Creates a dissertation-ready figure comparing VAF distributions:
  - Mixed strains (POS_MIXED_STRAINS)
  - Single-strain / other positives (other POS_* groups)

Adds a simple, clinically intuitive statistical support:
  Test 2: Enrichment of variants in an "intermediate AF band" (default 0.2–0.8)
    - Reports fraction in band for each bucket
    - Fisher's exact test (two-sided) on 2x2 table
    - Odds ratio + p-value

Inputs:
  --metadata  metadata.tsv (must contain 'sample' and 'group')
  --runs      one or more pipeline run output directories containing vaf_tables/
  --outdir    output directory

Outputs:
  Fig3_5_VAF_mixed_vs_single.png
  Table3_5_VAF_summary.tsv
  merged_vaf_for_3_5.tsv
"""

import argparse
from pathlib import Path
import sys
import re
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------
# Helpers
# --------------------------
def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of {candidates} found.\n"
        f"Available columns: {list(df.columns)}"
    )


def normalise_id(x: str) -> str:
    """
    Normalise sample IDs so that pipeline output filenames match metadata 'sample'.

    Examples:
      NEG_C013_sorted_vaf.tsv -> NEG_C013
      ROC001_sorted.tsv       -> ROC001
      MIX003_sorted_R1.tsv    -> MIX003
      DEF_POS_sorted.csv      -> DEF_POS
    """
    s = str(x)

    # drop extensions
    s = re.sub(r"\.(tsv|csv)$", "", s, flags=re.IGNORECASE)

    # remove common suffix words
    s = re.sub(r"(_vaf_table|_vaf)$", "", s, flags=re.IGNORECASE)

    # split at read marker if present
    s = re.split(r"_R1\b|_R2\b", s)[0]

    # split at _sorted (common pattern in your outputs)
    s = s.split("_sorted")[0]

    # split at _indelqual if present
    s = s.split("_indelqual")[0]

    # final cleanup
    s = s.strip("_- .")

    return s


def find_vaf_files(run_dir: Path):
    vdir = run_dir / "vaf_tables"
    if not vdir.is_dir():
        return []
    return sorted(list(vdir.glob("*.tsv")) + list(vdir.glob("*.csv")))


def read_vaf_file(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)

    # Try common AF/VAF column names
    af_col = None
    for cand in ["AF", "af", "VAF", "vaf", "allele_frequency", "alt_freq", "ALT_AF", "alt_af"]:
        if cand in df.columns:
            af_col = cand
            break
    if af_col is None:
        # last resort: pick something that looks like AF
        af_col = pick_col(df, ["ALT_FREQ", "alt_freq", "ALTFRAC", "altfrac"])

    out = pd.DataFrame()
    out["af"] = pd.to_numeric(df[af_col], errors="coerce")

    # base/sample ID from filename
    out["sample"] = normalise_id(path.name)
    out["_source_file"] = str(path)

    return out.dropna(subset=["af"])


def bucket_group(group: str) -> str:
    g = str(group)
    if g == "POS_MIXED_STRAINS":
        return "Mixed strains"
    if g.startswith("POS_"):
        return "Single strain / other positives"
    return "Ignore"


def frac_in_band(x: np.ndarray, lo: float, hi: float) -> float:
    if len(x) == 0:
        return float("nan")
    x = np.asarray(x, dtype=float)
    return float(((x >= lo) & (x <= hi)).mean())


# --------------------------
# Fisher exact (two-sided) WITHOUT SciPy
# --------------------------
def _log_choose(n: int, k: int) -> float:
    # log(C(n,k))
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_p(a: int, b: int, c: int, d: int) -> float:
    """
    Probability of observing table:
        [[a, b],
         [c, d]]
    given fixed margins, using hypergeometric.
    """
    r1 = a + b
    r2 = c + d
    c1 = a + c
    n = r1 + r2
    # p = [C(c1, a) * C(n - c1, r1 - a)] / C(n, r1)
    logp = _log_choose(c1, a) + _log_choose(n - c1, r1 - a) - _log_choose(n, r1)
    return float(math.exp(logp))


def fisher_exact_two_sided(a: int, b: int, c: int, d: int):
    """
    Returns (odds_ratio, p_value_two_sided) for a 2x2 table:
        [[a, b],
         [c, d]]
    Two-sided p is computed by summing probabilities of all tables
    with probability <= observed probability (classic Fisher two-sided).
    """
    # Odds ratio with a small guard
    odds_ratio = float("inf") if (b * c) == 0 and (a * d) > 0 else (a * d) / (b * c) if (b * c) > 0 else 0.0

    r1 = a + b
    r2 = c + d
    c1 = a + c
    c2 = b + d
    n = r1 + r2

    # Range of feasible a values given margins
    a_min = max(0, r1 - c2)
    a_max = min(r1, c1)

    p_obs = _hypergeom_p(a, b, c, d)

    p_two = 0.0
    for a_i in range(a_min, a_max + 1):
        b_i = r1 - a_i
        c_i = c1 - a_i
        d_i = r2 - c_i
        if min(b_i, c_i, d_i) < 0:
            continue
        p_i = _hypergeom_p(a_i, b_i, c_i, d_i)
        if p_i <= p_obs + 1e-15:
            p_two += p_i

    p_two = min(1.0, max(0.0, p_two))
    return odds_ratio, p_two


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot VAF distributions for mixed vs single-strain positives + intermediate-band Fisher test.")
    ap.add_argument("--metadata", required=True, help="metadata.tsv (truth) containing sample + group columns")
    ap.add_argument("--runs", required=True, nargs="+", help="Pipeline run output directories (contain vaf_tables/)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--dpi", type=int, default=450)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-points", type=int, default=25000)
    ap.add_argument("--no-points", action="store_true")

    ap.add_argument("--intermediate-lo", type=float, default=0.20)
    ap.add_argument("--intermediate-hi", type=float, default=0.80)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- load metadata
    meta = pd.read_csv(args.metadata, sep="\t")
    col_sample = pick_col(meta, ["sample", "base", "sample_id", "id"])
    col_group = pick_col(meta, ["group", "truth_group", "label_group"])
    meta = meta.rename(columns={col_sample: "sample", col_group: "group"})[["sample", "group"]].copy()
    meta["sample"] = meta["sample"].astype(str).apply(normalise_id)

    # --- load VAF tables from runs
    frames = []
    for r in args.runs:
        rdir = Path(r)
        if not rdir.is_dir():
            print(f"⚠️ Skipping (not a directory): {rdir}", file=sys.stderr)
            continue
        vaf_files = find_vaf_files(rdir)
        if not vaf_files:
            print(f"⚠️ No VAF files found in: {rdir}/vaf_tables", file=sys.stderr)
            continue
        for vf in vaf_files:
            try:
                frames.append(read_vaf_file(vf))
            except Exception as e:
                print(f"⚠️ Failed reading {vf}: {e}", file=sys.stderr)

    if not frames:
        sys.exit("❌ No VAF files loaded. Check your run directories contain vaf_tables/*.tsv or *.csv")

    vaf = pd.concat(frames, ignore_index=True)

    # --- join
    merged = vaf.merge(meta, on="sample", how="inner")

    if merged.empty:
        # Debug report to show what's not matching
        meta_ids = set(meta["sample"].astype(str))
        vaf_ids = set(vaf["sample"].astype(str))
        print("❌ Join produced 0 rows.", file=sys.stderr)
        print(f"Metadata sample IDs (example): {sorted(list(meta_ids))[:10]}", file=sys.stderr)
        print(f"VAF-derived sample IDs (example): {sorted(list(vaf_ids))[:10]}", file=sys.stderr)
        print(f"Metadata-only IDs (example): {sorted(list(meta_ids - vaf_ids))[:10]}", file=sys.stderr)
        print(f"VAF-only IDs (example): {sorted(list(vaf_ids - meta_ids))[:10]}", file=sys.stderr)
        sys.exit("Check naming: are your VAF filenames built from sample names or something else?")

    merged["bucket"] = merged["group"].apply(bucket_group)
    merged = merged.loc[merged["bucket"].isin(["Mixed strains", "Single strain / other positives"])].copy()

    if merged.empty:
        sys.exit("❌ Join succeeded but there are no POS_* samples with VAF tables in these runs.")

    merged_out = outdir / "merged_vaf_for_3_5.tsv"
    merged.to_csv(merged_out, sep="\t", index=False)

    # Downsample for plotting only (stats are done on full merged)
    merged_plot = merged
    if len(merged_plot) > args.max_points:
        merged_plot = merged_plot.sample(n=args.max_points, random_state=args.seed)

    # --------------------------
    # Test 2: intermediate-band enrichment + Fisher exact
    # --------------------------
    lo, hi = float(args.intermediate_lo), float(args.intermediate_hi)

    single_af = merged.loc[merged["bucket"] == "Single strain / other positives", "af"].astype(float).values
    mixed_af = merged.loc[merged["bucket"] == "Mixed strains", "af"].astype(float).values

    single_in = int(((single_af >= lo) & (single_af <= hi)).sum())
    single_out = int(((single_af < lo) | (single_af > hi)).sum())

    mixed_in = int(((mixed_af >= lo) & (mixed_af <= hi)).sum())
    mixed_out = int(((mixed_af < lo) | (mixed_af > hi)).sum())

    p_single = frac_in_band(single_af, lo, hi)
    p_mixed = frac_in_band(mixed_af, lo, hi)

    odds_ratio, p_fisher = fisher_exact_two_sided(mixed_in, mixed_out, single_in, single_out)

    # --------------------------
    # Summary table (includes Fisher test)
    # --------------------------
    rows = []
    for bucket in ["Single strain / other positives", "Mixed strains"]:
        sub = merged.loc[merged["bucket"] == bucket]
        af = sub["af"].astype(float).values
        n = len(af)
        pct_inter = float(((af >= lo) & (af <= hi)).mean() * 100.0) if n else 0.0
        rows.append({
            "bucket": bucket,
            "n_variants_total": int(n),
            f"percent_intermediate_AF_{lo:.2f}_to_{hi:.2f}": pct_inter,
            "median_AF": float(np.nanmedian(af)) if n else np.nan,
            "q25_AF": float(np.nanpercentile(af, 25)) if n else np.nan,
            "q75_AF": float(np.nanpercentile(af, 75)) if n else np.nan,
        })

    summary = pd.DataFrame(rows)

    # add a small "test" block as extra rows (easy to paste into Results)
    test_block = pd.DataFrame([{
        "bucket": "TEST_intermediate_band_fisher",
        "n_variants_total": int(len(single_af) + len(mixed_af)),
        f"percent_intermediate_AF_{lo:.2f}_to_{hi:.2f}": np.nan,
        "median_AF": np.nan,
        "q25_AF": np.nan,
        "q75_AF": np.nan,
        "single_in_band": single_in,
        "single_outside_band": single_out,
        "mixed_in_band": mixed_in,
        "mixed_outside_band": mixed_out,
        "p_single_in_band": float(p_single),
        "p_mixed_in_band": float(p_mixed),
        "odds_ratio_mixed_vs_single": float(odds_ratio),
        "fisher_p_two_sided": float(p_fisher),
        "band_lo": lo,
        "band_hi": hi,
    }])

    summary_out = pd.concat([summary, test_block], ignore_index=True)
    summary_path = outdir / "Table3_5_VAF_summary.tsv"
    summary_out.to_csv(summary_path, sep="\t", index=False)

    # --------------------------
    # Plot (cleaner, publication style)
    # --------------------------
    fig = plt.figure(figsize=(10.5, 6.0))
    ax = plt.gca()

    labels = ["Single strain / other positives", "Mixed strains"]
    data = [merged_plot.loc[merged_plot["bucket"] == lab, "af"].astype(float).values for lab in labels]

    # subtle intermediate band (keep shading, remove label text)
    ax.axhspan(lo, hi, alpha=0.06, zorder=0)

    vp = ax.violinplot(
        data,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        widths=0.82,
    )

    for body in vp["bodies"]:
        body.set_facecolor("#D9D9D9")
        body.set_edgecolor("#222222")
        body.set_linewidth(1.0)
        body.set_alpha(0.92)

    if not args.no_points:
        for i, lab in enumerate(labels, start=1):
            y = merged_plot.loc[merged_plot["bucket"] == lab, "af"].astype(float).values
            x = rng.normal(loc=i, scale=0.06, size=len(y))
            ax.scatter(
                x, y,
                s=8,
                alpha=0.18,
                color="black",
                linewidths=0,
                zorder=2,
            )

    ax.set_title("Variant allele frequency patterns distinguish mixed-strain from single-strain CMV samples", pad=10)
    ax.set_ylabel("Variant allele frequency (AF)")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Single strain / other positives", "Mixed strains"], rotation=12, ha="right")
    ax.set_ylim(0, 1.0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.15)

    # Optional: tiny caption line with stats (kept subtle)
    ax.text(
        0.02, 0.02,
        f"Intermediate band {lo:.2f}–{hi:.2f}: single={p_single:.3f}, mixed={p_mixed:.3f} | Fisher p={p_fisher:.2e}",
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=9,
        alpha=0.85,
    )

    plt.tight_layout()

    fig_path = outdir / "Fig3_5_VAF_mixed_vs_single.png"
    plt.savefig(fig_path, dpi=args.dpi)
    plt.close(fig)

    # --------------------------
    # Print a clean Results-friendly summary
    # --------------------------
    print("✅ Done")
    print(f"Saved figure: {fig_path}")
    print(f"Saved table:  {summary_path}")
    print(f"Saved merged: {merged_out}\n")

    print("Intermediate AF band enrichment (Test 2)")
    print(f"Band: {lo:.2f}–{hi:.2f}")
    print(f"Single strain / other positives: {single_in}/{single_in + single_out} = {p_single:.3f}")
    print(f"Mixed strains:                  {mixed_in}/{mixed_in + mixed_out} = {p_mixed:.3f}")
    print(f"Fisher's exact (two-sided): OR={odds_ratio:.3f}, p={p_fisher:.3e}")


if __name__ == "__main__":
    main()