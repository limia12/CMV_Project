#!/usr/bin/env python3
import argparse
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------
# Naming / canonicalization
# ----------------------------
_SUFFIX_PATTERNS = [
    r'(?:_R?[12])$',
    r'(?:_00[1-9])$',
    r'(?:_sorted(?:_sorted)*)$',
    r'(?:_sort(?:ed)?)$',
    r'(?:_dedup(?:ed)?)$',
    r'(?:_mark(?:ed)?dups?)$',
    r'(?:_rmdup[s]?)$',
    r'(?:_trim(?:med)?(?:_paired)?)$',
    r'(?:_filtered)$',
    r'(?:\.bam)$',
]

def canonical_sample(name: str) -> str:
    s = name
    changed = True
    while changed:
        changed = False
        for pat in _SUFFIX_PATTERNS:
            new = re.sub(pat, "", s, flags=re.IGNORECASE)
            if new != s:
                s = new
                changed = True
    return s.rstrip("_-")

def base_from_blast_filename(p: Path) -> str:
    # e.g. COINF01_sorted_R1_blast.tsv -> COINF01
    n = p.name
    n = re.sub(r'_blast\.tsv$', '', n, flags=re.IGNORECASE)
    n = re.sub(r'_(R1|R2)$', '', n, flags=re.IGNORECASE)
    return canonical_sample(n)

def base_from_vaf_filename(p: Path) -> str:
    # e.g. COINF01_sorted_vaf.tsv -> COINF01
    n = p.name
    n = re.sub(r'_vaf\.tsv$', '', n, flags=re.IGNORECASE)
    return canonical_sample(n)


# ----------------------------
# BLAST parsing (streaming)
# ----------------------------
CMV_PAT = re.compile(
    r"Human\s+herpesvirus\s*5|Human\s+cytomegalovirus|Cytomegalovirus|HCMV|NC_006273\.2|Merlin",
    re.IGNORECASE
)

def safe_evalue(x: float) -> float:
    if not math.isfinite(x) or x <= 0:
        return 1e-300
    return x

def parse_blast_tsv(path: Path) -> Dict[str, float]:
    """
    Reads outfmt:
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    Returns totals and CMV-only summary stats.
    """
    total = 0
    cmv = 0
    non_cmv = 0

    sum_pident = 0.0
    cmv_pident_n = 0

    # for geometric mean: log(evalue)
    sum_log_e = 0.0
    cmv_e_n = 0

    # mean -log10(e)
    sum_neglog10 = 0.0
    cmv_nl10_n = 0

    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip():
                continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            # guard: some lines may be malformed
            if len(parts) < 13:
                continue

            try:
                pident = float(parts[2])
            except Exception:
                pident = float("nan")

            try:
                evalue = float(parts[10])
            except Exception:
                evalue = float("nan")

            stitle = parts[12] if len(parts) > 12 else ""

            is_cmv = bool(CMV_PAT.search(stitle))
            if is_cmv:
                cmv += 1
                if math.isfinite(pident):
                    sum_pident += pident
                    cmv_pident_n += 1

                ev = safe_evalue(evalue)
                sum_log_e += math.log(ev)
                cmv_e_n += 1

                sum_neglog10 += (-math.log10(ev))
                cmv_nl10_n += 1
            else:
                non_cmv += 1

    mean_pident = (sum_pident / cmv_pident_n) if cmv_pident_n else float("nan")
    geo_mean_e = math.exp(sum_log_e / cmv_e_n) if cmv_e_n else float("nan")
    mean_neglog10 = (sum_neglog10 / cmv_nl10_n) if cmv_nl10_n else float("nan")

    return {
        "blast_total_hits": float(total),
        "blast_cmv_hits": float(cmv),
        "blast_non_cmv_hits": float(non_cmv),
        "blast_mean_cmv_pident": mean_pident,
        "blast_geo_mean_cmv_evalue": geo_mean_e,
        "blast_mean_neglog10_cmv_evalue": mean_neglog10,
    }

def aggregate_blast(blast_dir: Path) -> pd.DataFrame:
    files = sorted(blast_dir.glob("*_blast.tsv"))
    if not files:
        return pd.DataFrame(columns=[
            "sample","blast_total_hits","blast_cmv_hits","blast_non_cmv_hits",
            "blast_mean_cmv_pident","blast_geo_mean_cmv_evalue","blast_mean_neglog10_cmv_evalue"
        ])

    # sum counts across R1+R2 per base, and combine means weighted by CMV hit counts
    accum: Dict[str, Dict[str, float]] = {}

    # For weighted means we need numerator+denominator; easiest is re-parse with cmv count,
    # but we didn't return cmv count separately. We'll approximate using blast_cmv_hits as weight.
    for f in files:
        base = base_from_blast_filename(f)
        stats = parse_blast_tsv(f)
        if base not in accum:
            accum[base] = {
                "blast_total_hits": 0.0,
                "blast_cmv_hits": 0.0,
                "blast_non_cmv_hits": 0.0,
                "pident_sum_w": 0.0,
                "pident_w": 0.0,
                "neglog10_sum_w": 0.0,
                "neglog10_w": 0.0,
                "log_e_sum_w": 0.0,
                "log_e_w": 0.0,
            }

        a = accum[base]
        a["blast_total_hits"] += stats["blast_total_hits"]
        a["blast_cmv_hits"] += stats["blast_cmv_hits"]
        a["blast_non_cmv_hits"] += stats["blast_non_cmv_hits"]

        w = stats["blast_cmv_hits"]
        if math.isfinite(stats["blast_mean_cmv_pident"]) and w > 0:
            a["pident_sum_w"] += stats["blast_mean_cmv_pident"] * w
            a["pident_w"] += w

        if math.isfinite(stats["blast_mean_neglog10_cmv_evalue"]) and w > 0:
            a["neglog10_sum_w"] += stats["blast_mean_neglog10_cmv_evalue"] * w
            a["neglog10_w"] += w

        if math.isfinite(stats["blast_geo_mean_cmv_evalue"]) and w > 0:
            # geo mean combine via logs:
            a["log_e_sum_w"] += math.log(safe_evalue(stats["blast_geo_mean_cmv_evalue"])) * w
            a["log_e_w"] += w

    rows = []
    for base, a in accum.items():
        mean_pident = a["pident_sum_w"] / a["pident_w"] if a["pident_w"] else float("nan")
        mean_neglog10 = a["neglog10_sum_w"] / a["neglog10_w"] if a["neglog10_w"] else float("nan")
        geo_mean_e = math.exp(a["log_e_sum_w"] / a["log_e_w"]) if a["log_e_w"] else float("nan")

        rows.append({
            "sample": base,
            "blast_total_hits": int(a["blast_total_hits"]),
            "blast_cmv_hits": int(a["blast_cmv_hits"]),
            "blast_non_cmv_hits": int(a["blast_non_cmv_hits"]),
            "blast_mean_cmv_pident": mean_pident,
            "blast_geo_mean_cmv_evalue": geo_mean_e,
            "blast_mean_neglog10_cmv_evalue": mean_neglog10,
        })
    return pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)


# ----------------------------
# VAF parsing
# ----------------------------
def aggregate_vaf(vaf_dir: Path, mid_lo: float, mid_hi: float) -> pd.DataFrame:
    files = sorted(vaf_dir.glob("*_vaf.tsv"))
    if not files:
        return pd.DataFrame(columns=["sample","vaf_total_variants","vaf_mid_af_count","vaf_mid_af_frac"])

    rows = []
    for f in files:
        base = base_from_vaf_filename(f)
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue

        # Expect columns: CHROM POS REF ALT DP AD_ALT AF
        if "AF" not in df.columns:
            # sometimes AF may be string; try coerce
            pass

        af = pd.to_numeric(df.get("AF", pd.Series([], dtype=float)), errors="coerce").dropna()
        total = int(len(af))
        mid = int(((af >= mid_lo) & (af <= mid_hi)).sum()) if total else 0
        frac = float(mid / total) if total else float("nan")

        rows.append({
            "sample": base,
            "vaf_total_variants": total,
            "vaf_mid_af_count": mid,
            "vaf_mid_af_frac": frac,
        })

    return pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)


# ----------------------------
# Mutation detection vs truth
# ----------------------------
_MUT_APPLIED_RE = re.compile(r"^(\d+)([ACGTN])>([ACGTN])$", re.IGNORECASE)

def parse_truth_mutations(muts_applied: str) -> List[Tuple[int, str, str]]:
    """
    metadata.tsv column "mutations_applied" contains comma-separated tokens like:
      141800A>G,78200C>T
    Returns list of (pos, ref, alt)
    """
    out = []
    if not isinstance(muts_applied, str) or not muts_applied.strip():
        return out
    for tok in muts_applied.split(","):
        tok = tok.strip()
        if not tok:
            continue
        m = _MUT_APPLIED_RE.match(tok)
        if not m:
            continue
        pos = int(m.group(1))
        ref = m.group(2).upper()
        alt = m.group(3).upper()
        out.append((pos, ref, alt))
    return out

def detect_mutations_from_vaf(truth_meta: pd.DataFrame, vaf_dir: Path, min_detect_af: float) -> pd.DataFrame:
    """
    For each mutation-spiked sample, check whether truth mutations appear in its *_vaf.tsv.
    """
    # load all vaf tables into dict for quick lookup
    vafs: Dict[str, pd.DataFrame] = {}
    for f in sorted(vaf_dir.glob("*_vaf.tsv")):
        base = base_from_vaf_filename(f)
        try:
            df = pd.read_csv(f, sep="\t")
            df["POS"] = pd.to_numeric(df.get("POS"), errors="coerce")
            df["AF"] = pd.to_numeric(df.get("AF"), errors="coerce")
            df["ALT"] = df.get("ALT", "").astype(str)
            vafs[base] = df
        except Exception:
            continue

    rows = []
    for _, r in truth_meta.iterrows():
        sample = str(r["sample"])
        group = str(r.get("group", ""))
        mutation_spiked = int(r.get("mutation_spiked", 0)) == 1
        target_af = r.get("mutation_af_target", "")
        muts = parse_truth_mutations(str(r.get("mutations_applied", "")))
        if not mutation_spiked:
            continue

        df = vafs.get(sample, pd.DataFrame())
        for (pos, ref, alt) in muts:
            found = 0
            best_af = np.nan
            if not df.empty:
                sub = df[(df["POS"] == pos) & (df["ALT"].astype(str) == alt)].copy()
                if not sub.empty:
                    best_af = float(np.nanmax(sub["AF"].to_numpy(dtype=float)))
                    if np.isfinite(best_af) and best_af >= min_detect_af:
                        found = 1

            rows.append({
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "found": found,
                "af": best_af if np.isfinite(best_af) else "",
                "sample": sample,
                "group": group,
                "target_af": target_af,
                "detected_af_ge_min": 1 if found == 1 else 0
            })

    if not rows:
        return pd.DataFrame(columns=["pos","ref","alt","found","af","sample","group","target_af","detected_af_ge_min"])
    return pd.DataFrame(rows).sort_values(["sample","pos"]).reset_index(drop=True)


# ----------------------------
# Simple statistics (self-contained)
# ----------------------------
def confusion_counts(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[int,int,int,int]:
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn

def safe_div(a: float, b: float) -> float:
    return float(a / b) if b else float("nan")

def fisher_exact_2x2(tp: int, fp: int, fn: int, tn: int) -> float:
    """
    Two-sided Fisher exact p-value via hypergeometric enumeration.
    Table layout:
        pred+  pred-
    true+   tp    fn
    true-   fp    tn
    """
    # Margins:
    r1 = tp + fn
    r2 = fp + tn
    c1 = tp + fp
    c2 = fn + tn
    n = r1 + r2

    def hypergeom_p(x: int) -> float:
        # probability of x in (true+, pred+) given margins
        return (math.comb(r1, x) * math.comb(r2, c1 - x)) / math.comb(n, c1)

    # feasible x range
    lo = max(0, c1 - r2)
    hi = min(r1, c1)
    p_obs = hypergeom_p(tp)
    p = 0.0
    for x in range(lo, hi + 1):
        px = hypergeom_p(x)
        if px <= p_obs + 1e-12:
            p += px
    return min(1.0, p)

def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """
    Spearman rho via rank correlation (average ranks for ties).
    """
    def rankdata(a: np.ndarray) -> np.ndarray:
        order = np.argsort(a, kind="mergesort")
        ranks = np.empty(len(a), dtype=float)
        ranks[order] = np.arange(1, len(a) + 1, dtype=float)
        # average ties
        _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
        for k, c in enumerate(counts):
            if c > 1:
                idx = np.where(inv == k)[0]
                ranks[idx] = ranks[idx].mean()
        return ranks

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]; y = y[mask]
    if len(x) < 3:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    # Pearson on ranks
    rxm = rx - rx.mean()
    rym = ry - ry.mean()
    denom = math.sqrt((rxm**2).sum() * (rym**2).sum())
    return float((rxm * rym).sum() / denom) if denom else float("nan")


# ----------------------------
# Plot helpers (matplotlib only)
# ----------------------------
def save_confusion_matrix(tp, tn, fp, fn, out_png: Path):
    mat = np.array([[tp, fn],
                    [fp, tn]], dtype=int)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.imshow(mat)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["Pred +", "Pred -"])
    ax.set_yticklabels(["True +", "True -"])
    ax.set_title("CMV detection confusion matrix")

    for (i, j), v in np.ndenumerate(mat):
        ax.text(j, i, str(v), ha="center", va="center")

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def save_boxplot_by_group(df: pd.DataFrame, value_col: str, group_col: str, title: str, ylabel: str, out_png: Path):
    g = df[[group_col, value_col]].copy()
    g = g[g[value_col].notna()]
    groups = list(g[group_col].unique())
    data = [g[g[group_col] == k][value_col].to_numpy(dtype=float) for k in groups]

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.boxplot(data, labels=groups, showfliers=True)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def save_scatter(df: pd.DataFrame, x: str, y: str, title: str, xlabel: str, ylabel: str, out_png: Path):
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.scatter(df[x].to_numpy(dtype=float), df[y].to_numpy(dtype=float))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def save_threshold_curve(y_true: np.ndarray, scores: np.ndarray, out_png: Path):
    """
    Plot sensitivity & specificity across a range of BLAST CMV-hit thresholds.
    """
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    thresholds = np.unique(np.clip(scores[~np.isnan(scores)], 0, None).astype(int))
    if thresholds.size == 0:
        thresholds = np.arange(0, 101, 5)
    else:
        thresholds = np.unique(np.concatenate([thresholds, np.arange(0, int(thresholds.max()) + 1, max(1, int(thresholds.max() / 20)))]))

    sens = []
    spec = []
    for t in thresholds:
        y_pred = (scores >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y_true, y_pred)
        sens.append(safe_div(tp, tp + fn))
        spec.append(safe_div(tn, tn + fp))

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.plot(thresholds, sens, label="Sensitivity (Recall)")
    ax.plot(thresholds, spec, label="Specificity")
    ax.set_title("Detection performance vs threshold (BLAST CMV-hit cutoff)")
    ax.set_xlabel("Threshold: CMV hits ≥ t")
    ax.set_ylabel("Rate")
    ax.set_ylim(-0.02, 1.02)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def save_mutation_target_vs_observed(mut_df: pd.DataFrame, out_png: Path):
    if mut_df.empty:
        return
    d = mut_df.copy()
    d["target_af"] = pd.to_numeric(d["target_af"], errors="coerce")
    d["af_num"] = pd.to_numeric(d["af"], errors="coerce")

    # per sample: best observed AF across truth mutations
    per = d.groupby(["sample","target_af"], as_index=False)["af_num"].max()
    per = per.dropna(subset=["target_af"])

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.scatter(per["target_af"].to_numpy(dtype=float), per["af_num"].to_numpy(dtype=float))
    ax.set_title("Mutation spike: target AF vs observed max AF (from VAF tables)")
    ax.set_xlabel("Target AF (truth)")
    ax.set_ylabel("Observed max AF (pipeline)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(description="Evaluate CMV pipeline on synthetic stress-test set and generate figures.")
    ap.add_argument("--truth", required=True, help="Path to metadata.tsv produced by the ART generator.")
    ap.add_argument("--blast-dir", required=True, help="Folder containing *_blast.tsv (R1/R2) outputs.")
    ap.add_argument("--vaf-dir", required=True, help="Folder containing *_vaf.tsv outputs.")
    ap.add_argument("--outdir", required=True, help="Output directory for tables + figures.")
    ap.add_argument("--min-cmv-hits", type=int, default=25, help="Threshold for calling CMV detected from BLAST CMV hits.")
    ap.add_argument("--mid-af-lo", type=float, default=0.2, help="Lower AF bound for 'mid-AF' mixed-strain heuristic.")
    ap.add_argument("--mid-af-hi", type=float, default=0.8, help="Upper AF bound for 'mid-AF' mixed-strain heuristic.")
    ap.add_argument("--min-detect-af", type=float, default=0.01, help="Minimum AF to count a truth mutation as detected.")
    args = ap.parse_args()

    truth_path = Path(args.truth)
    blast_dir = Path(args.blast_dir)
    vaf_dir = Path(args.vaf_dir)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    # ---- Load truth ----
    truth = pd.read_csv(truth_path, sep="\t", dtype=str).fillna("")
    if "sample" not in truth.columns:
        raise SystemExit(f"truth file missing 'sample' column: {truth_path}")

    # Ensure expected columns exist (some are optional depending on your generator version)
    for c in ["group", "cmv_present", "mutation_spiked", "mutation_af_target", "mutations_applied"]:
        if c not in truth.columns:
            truth[c] = ""

    # Normalize types
    truth["sample"] = truth["sample"].astype(str).apply(canonical_sample)
    truth["truth_cmv_present"] = truth["cmv_present"].astype(str).replace({"": "0"}).astype(int)
    truth["mutation_spiked"] = truth["mutation_spiked"].astype(str).replace({"": "0"}).astype(int)

    # expected fraction: from generator metadata it may be embedded differently; we keep as optional
    if "expected_cmv_fraction" not in truth.columns:
        # if you have cmv_pairs/human_pairs etc you could compute it — left as blank unless present
        truth["expected_cmv_fraction"] = ""

    # ---- Aggregate pipeline outputs ----
    blast_df = aggregate_blast(blast_dir)
    vaf_df = aggregate_vaf(vaf_dir, mid_lo=args.mid_af_lo, mid_hi=args.mid_af_hi)

    # ---- Join ----
    df = truth.merge(blast_df, on="sample", how="left").merge(vaf_df, on="sample", how="left")
    for c in ["blast_total_hits","blast_cmv_hits","blast_non_cmv_hits","vaf_total_variants","vaf_mid_af_count"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    if "vaf_mid_af_frac" in df.columns:
        df["vaf_mid_af_frac"] = pd.to_numeric(df["vaf_mid_af_frac"], errors="coerce")

    # ---- CMV detection call ----
    df["pred_cmv_detected"] = (df["blast_cmv_hits"] >= int(args.min_cmv_hits)).astype(int)

    # ---- Mutation detection table ----
    mut_df = detect_mutations_from_vaf(truth, vaf_dir, min_detect_af=float(args.min_detect_af))
    mut_df.to_csv(outdir / "mutation_table.tsv", sep="\t", index=False)

    # ---- Save evaluation table ----
    # Keep a tidy subset + keep everything else if you want to expand
    keep_cols = [
        "sample","group","truth_cmv_present","pred_cmv_detected","expected_cmv_fraction",
        "blast_total_hits","blast_cmv_hits","blast_non_cmv_hits",
        "blast_mean_cmv_pident","blast_geo_mean_cmv_evalue","blast_mean_neglog10_cmv_evalue",
        "vaf_total_variants","vaf_mid_af_count","vaf_mid_af_frac",
        "mutation_spiked","mutation_af_target"
    ]
    for c in keep_cols:
        if c not in df.columns:
            df[c] = ""
    df[keep_cols].to_csv(outdir / "evaluation_table.tsv", sep="\t", index=False)

    # ---- Stats: detection performance ----
    y_true = df["truth_cmv_present"].to_numpy(dtype=int)
    y_pred = df["pred_cmv_detected"].to_numpy(dtype=int)
    tp, tn, fp, fn = confusion_counts(y_true, y_pred)

    sens = safe_div(tp, tp + fn)
    spec = safe_div(tn, tn + fp)
    prec = safe_div(tp, tp + fp)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    acc = safe_div(tp + tn, tp + tn + fp + fn)

    p_fisher = fisher_exact_2x2(tp, fp, fn, tn)

    # Spearman: expected fraction vs CMV hits (only if expected fraction present)
    frac = pd.to_numeric(df["expected_cmv_fraction"], errors="coerce")
    rho = spearman_rho(frac.to_numpy(dtype=float), df["blast_cmv_hits"].to_numpy(dtype=float))

    # Mutation detection rate summary
    mut_rate = float(mut_df["found"].mean()) if not mut_df.empty else float("nan")

    # ---- Write stats summary ----
    with open(outdir / "stats_summary.txt", "w") as fh:
        fh.write("CMV detection evaluation (truth vs pipeline)\n")
        fh.write(f"min_cmv_hits threshold: {args.min_cmv_hits}\n\n")
        fh.write(f"TP={tp}  TN={tn}  FP={fp}  FN={fn}\n")
        fh.write(f"Sensitivity/Recall: {sens:.3f}\n")
        fh.write(f"Specificity:         {spec:.3f}\n")
        fh.write(f"Precision:           {prec:.3f}\n")
        fh.write(f"F1:                  {f1:.3f}\n")
        fh.write(f"Accuracy:            {acc:.3f}\n\n")
        fh.write("Significance tests\n")
        fh.write(f"Fisher exact (2-sided) p-value for association (truth vs prediction): {p_fisher:.6g}\n")
        fh.write(f"Spearman rho (expected CMV fraction vs BLAST CMV hits): {rho:.3f}\n\n")
        fh.write("Mutation detection\n")
        fh.write(f"Per-truth-mutation detection rate (AF >= {args.min_detect_af}): {mut_rate:.3f}\n")

    # ---- Figures ----
    save_confusion_matrix(tp, tn, fp, fn, figdir / "fig01_confusion_matrix.png")

    save_threshold_curve(
        y_true=y_true,
        scores=df["blast_cmv_hits"].to_numpy(dtype=float),
        out_png=figdir / "fig02_threshold_curve.png"
    )

    # Boxplot of CMV hits by group
    if "group" in df.columns and df["group"].astype(str).str.len().sum() > 0:
        save_boxplot_by_group(
            df=df,
            value_col="blast_cmv_hits",
            group_col="group",
            title="BLAST CMV hit counts by truth group",
            ylabel="CMV hits (BLAST lines matching CMV title pattern)",
            out_png=figdir / "fig03_blast_cmv_hits_by_group.png"
        )

    # Mixed strain heuristic figure
    if "vaf_mid_af_frac" in df.columns and df["vaf_mid_af_frac"].notna().any():
        save_boxplot_by_group(
            df=df,
            value_col="vaf_mid_af_frac",
            group_col="group",
            title=f"Mixed-strain signal: fraction of mid-AF variants ({args.mid_af_lo}-{args.mid_af_hi}) by group",
            ylabel="Mid-AF fraction",
            out_png=figdir / "fig04_mid_af_fraction_by_group.png"
        )

    # Expected fraction vs CMV hits (only if fraction column has values)
    frac_num = pd.to_numeric(df["expected_cmv_fraction"], errors="coerce")
    tmp = df.copy()
    tmp["expected_cmv_fraction_num"] = frac_num
    tmp = tmp.dropna(subset=["expected_cmv_fraction_num"])
    if not tmp.empty:
        save_scatter(
            df=tmp,
            x="expected_cmv_fraction_num",
            y="blast_cmv_hits",
            title="Expected CMV fraction (truth) vs BLAST CMV hits (pipeline)",
            xlabel="Expected CMV fraction (truth)",
            ylabel="BLAST CMV hits",
            out_png=figdir / "fig05_expected_fraction_vs_hits.png"
        )

    # Mutation AF target vs observed
    if not mut_df.empty:
        save_mutation_target_vs_observed(mut_df, figdir / "fig06_mutation_target_vs_observed.png")

    print("\n✅ Wrote:")
    print(f"  - {outdir / 'evaluation_table.tsv'}")
    print(f"  - {outdir / 'mutation_table.tsv'}")
    print(f"  - {outdir / 'stats_summary.txt'}")
    print(f"  - figures in {figdir}")

if __name__ == "__main__":
    main()
