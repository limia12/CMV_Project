#!/usr/bin/env python3
"""
cmv_threshold_optimize.py

Purpose
-------
Pick a robust CMV detection threshold that uses *multiple* metrics (not just mapped HQ reads),
and generate dissertation-ready, well-annotated plots explaining why the chosen threshold works.

This script is designed to plug directly into your existing pipeline outputs:

1) Alignment summary (from cmv_pipeline_gui.py samtools metrics):
   - cmv_alignment_summary.csv (preferred) or *_summary.csv
   - columns include:
       base, cmv_mapped_reads_mapq, cmv_rpm_mapq, cmv_breadth_mapq,
       cmv_cov_bases_mapq, cmv_ref_bases, total_reads_in_cmv_bam

2) Gene-level coverage (from extract_cmv_bam_metrics.py):
   - *_gene.csv with columns including: base, gene, avg_depth, breadth

3) BLAST outputs (from summarize_blast_hits.py and/or raw *_blast.tsv):
   - *_blast_top5.csv with columns: base, sample, species, total_hits, avg_identity
   - optionally raw *_blast.tsv loaded (outfmt 6), columns include pident, stitle, length, bitscore

Truth
-----
Use the metadata.tsv produced by your simulator script (truth/metadata.tsv).
It contains:
  sample, group, truth_cmv_present, cmv_fraction_target, ...

Multi-directory support
-----------------------
If you had to run the pipeline in batches (e.g., 256 samples crashing), you can pass multiple
directories for CSV outputs (and optional BLAST hit dirs). The script will union them and
deduplicate on canonical 'base' sample names.

What it does (high-level)
-------------------------
A) Builds a modelling table at "base" sample level:
   - y_true (truth from metadata.tsv)
   - candidate metrics (mapped HQ reads, breadth, RPM, gene depth, BLAST identity, etc.)
B) Exploratory plots:
   - distribution by truth class, correlation heatmap, and per-metric ROC/PR curves
C) Statistical support:
   - permutation test for AUC (per metric and selected combined model)
   - bootstrap confidence intervals for AUC, AP, F1, recall, specificity, and precision
   - Bayesian posterior distributions for sensitivity/specificity/PPV/NPV at chosen threshold
D) Threshold selection:
   - compares *single-metric* thresholds vs *multi-metric* rules
   - selects the best approach by:
       1) maximize F1
       2) tie-break: higher recall
       3) tie-break: higher specificity
       4) tie-break: simpler model (fewer metrics)
E) Saves:
   - results/thresholds.json (final rule + numbers)
   - results/metrics_table.tsv (modelling table)
   - results/per_metric_performance.tsv
   - results/model_comparison.tsv
   - results/figures/*.png (publication-ready)

No sklearn dependency
---------------------
This script avoids scikit-learn to keep it lightweight and consistent with your ROC code style.
(Uses numpy/pandas/matplotlib only.)
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# canonical sample naming (match your pipeline style)
# -----------------------------
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
    """
    Convert a filename-ish sample string into a stable "base" sample name.
    Repeatedly strips common suffix patterns until nothing else can be removed.
    """
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


# -----------------------------
# IO helpers
# -----------------------------
def die(msg: str) -> None:
    print(f"❌ {msg}", file=sys.stderr)
    raise SystemExit(2)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_table(path: Path) -> pd.DataFrame:
    """
    Robust CSV/TSV reader with auto separator and safe dtype handling.
    """
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    for sep in [",", "\t"]:
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


def read_glob(paths: Sequence[Path]) -> pd.DataFrame:
    dfs = []
    for p in paths:
        df = read_table(p)
        if df is not None and not df.empty:
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def coerce_numeric(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# -----------------------------
# metric extraction from your outputs
# -----------------------------
def load_alignment_summaries(csv_dirs: List[Path]) -> pd.DataFrame:
    """
    Load cmv_alignment_summary.csv if present, otherwise *_summary.csv.
    Expects the samtools-derived columns used in cmv_pipeline_gui.py.
    """
    files = []
    for d in csv_dirs:
        if not d.exists():
            continue
        p = d / "cmv_alignment_summary.csv"
        if p.exists():
            files.append(p)
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
            guess = df.columns[0]
            df["base"] = df[guess].astype(str).apply(canonical_sample)

    df["base"] = df["base"].astype(str).apply(canonical_sample)

    numeric_cols = [
        "cmv_mapped_reads_mapq",
        "cmv_rpm_mapq",
        "cmv_breadth_mapq",
        "cmv_cov_bases_mapq",
        "cmv_ref_bases",
        "total_reads_in_cmv_bam",
    ]
    df = coerce_numeric(df, numeric_cols)

    agg = {"cmv_mapped_reads_mapq": "sum",
           "cmv_cov_bases_mapq": "sum",
           "total_reads_in_cmv_bam": "sum",
           "cmv_rpm_mapq": "max",
           "cmv_breadth_mapq": "max",
           "cmv_ref_bases": "max"}
    for k in list(agg.keys()):
        if k not in df.columns:
            agg.pop(k, None)

    out = df.groupby("base", as_index=False).agg(agg) if agg else df[["base"]].drop_duplicates()
    return out


def load_gene_metrics(csv_dirs: List[Path]) -> pd.DataFrame:
    """
    Load *_gene.csv outputs and summarise to base-level metrics.
    """
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
            die("Gene table missing 'base' and 'sample' columns.")

    df["base"] = df["base"].astype(str).apply(canonical_sample)
    df = coerce_numeric(df, ["avg_depth", "breadth"])

    g = df.groupby("base")
    base_df = pd.DataFrame({
        "base": list(g.groups.keys()),
        "gene_mean_depth": g["avg_depth"].mean().values,
        "gene_median_depth": g["avg_depth"].median().values,
        "gene_min_depth": g["avg_depth"].min().values,
        "gene_mean_breadth": g["breadth"].mean().values,
        "gene_min_breadth": g["breadth"].min().values,
        "n_genes_reported": g.size().values,
    })

    def gene_val(gene: str, col: str) -> pd.Series:
        sub = df[df["gene"].astype(str).str.upper() == gene.upper()]
        if sub.empty:
            return pd.Series(index=base_df["base"], data=np.nan)
        return sub.groupby("base")[col].mean()

    for gene in ["UL97", "UL54", "IE1", "UL83"]:
        s = gene_val(gene, "avg_depth")
        base_df[f"{gene}_mean_depth"] = base_df["base"].map(s).astype(float)

    return base_df


def load_blast_top5(csv_dirs: List[Path]) -> pd.DataFrame:
    """
    Load *_blast_top5.csv outputs and compute CMV-focused identity stats per base.
    """
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
            die("BLAST top5 table missing 'base' and 'sample' columns.")

    df["base"] = df["base"].astype(str).apply(canonical_sample)
    df = coerce_numeric(df, ["total_hits", "avg_identity"])

    species = df["species"].astype(str) if "species" in df.columns else pd.Series("", index=df.index)
    is_cmv = species.str.contains("herpesvirus 5|cytomegalovirus|merlin|hcmv", case=False, na=False)

    out = []
    for b, sub in df.groupby("base"):
        cmv_sub = sub[is_cmv.loc[sub.index]]
        if cmv_sub.empty:
            out.append({"base": b, "blast_cmv_mean_identity_top5": np.nan,
                        "blast_cmv_hits_top5": 0,
                        "blast_any_hits_top5": int(sub["total_hits"].sum()) if "total_hits" in sub else np.nan})
        else:
            if "total_hits" in cmv_sub.columns and cmv_sub["total_hits"].sum() > 0:
                w = cmv_sub["total_hits"].values
                mu = float(np.average(cmv_sub["avg_identity"].values, weights=w))
                hits = int(cmv_sub["total_hits"].sum())
            else:
                mu = float(cmv_sub["avg_identity"].mean())
                hits = int(cmv_sub.shape[0])
            out.append({"base": b, "blast_cmv_mean_identity_top5": mu,
                        "blast_cmv_hits_top5": hits,
                        "blast_any_hits_top5": int(sub["total_hits"].sum()) if "total_hits" in sub else np.nan})
    return pd.DataFrame(out)


def load_truth(metadata_tsv: Path) -> pd.DataFrame:
    df = read_table(metadata_tsv)
    if df.empty:
        die(f"Truth metadata empty or unreadable: {metadata_tsv}")
    required = ["sample", "truth_cmv_present"]
    for r in required:
        if r not in df.columns:
            die(f"Truth metadata missing required column '{r}'")
    df["base"] = df["sample"].astype(str).apply(canonical_sample)
    df["y_true"] = pd.to_numeric(df["truth_cmv_present"], errors="coerce").fillna(0).astype(int)
    keep = ["base", "y_true"]
    for c in ["group", "cmv_fraction_target", "total_pairs"]:
        if c in df.columns:
            keep.append(c)
    out = df[keep].drop_duplicates("base")
    return out


# -----------------------------
# metrics + scoring
# -----------------------------
@dataclass
class Perf:
    threshold: float
    tp: int
    tn: int
    fp: int
    fn: int
    precision: float
    recall: float
    specificity: float
    f1: float


def confusion_from_pred(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[int, int, int, int]:
    y_true = y_true.astype(int)
    y_pred = y_pred.astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn


def safe_div(a: float, b: float) -> float:
    return float(a / b) if b else 0.0


def perf_from_threshold(y_true: np.ndarray, x: np.ndarray, thr: float, direction: str = "ge") -> Perf:
    if direction == "ge":
        y_pred = (x >= thr).astype(int)
    else:
        y_pred = (x <= thr).astype(int)
    tp, tn, fp, fn = confusion_from_pred(y_true, y_pred)
    prec = safe_div(tp, tp + fp)
    rec = safe_div(tp, tp + fn)
    spec = safe_div(tn, tn + fp)
    f1 = safe_div(2 * prec * rec, prec + rec) if (prec + rec) else 0.0
    return Perf(threshold=float(thr), tp=tp, tn=tn, fp=fp, fn=fn,
                precision=prec, recall=rec, specificity=spec, f1=f1)


def roc_curve_points(y_true: np.ndarray, scores: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_true = y_true.astype(int)
    m = np.isfinite(scores)
    y = y_true[m]
    s = scores[m]
    thrs = np.unique(s)[::-1]
    thrs = np.r_[np.inf, thrs, -np.inf]
    tpr = []
    fpr = []
    for thr in thrs:
        y_pred = (s >= thr).astype(int)
        tp, tn, fp, fn = confusion_from_pred(y, y_pred)
        tpr.append(safe_div(tp, tp + fn))
        fpr.append(safe_div(fp, fp + tn))
    return np.array(fpr), np.array(tpr), thrs


def pr_curve_points(y_true: np.ndarray, scores: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_true = y_true.astype(int)
    m = np.isfinite(scores)
    y = y_true[m]
    s = scores[m]
    thrs = np.unique(s)[::-1]
    thrs = np.r_[np.inf, thrs, -np.inf]
    precs = []
    recs = []
    for thr in thrs:
        y_pred = (s >= thr).astype(int)
        tp, tn, fp, fn = confusion_from_pred(y, y_pred)
        precs.append(safe_div(tp, tp + fp))
        recs.append(safe_div(tp, tp + fn))
    return np.array(recs), np.array(precs), thrs


def auc_trapz(x: np.ndarray, y: np.ndarray) -> float:
    idx = np.argsort(x)
    xs = x[idx]
    ys = y[idx]
    return float(np.trapezoid(ys, xs))


def average_precision(y_true: np.ndarray, scores: np.ndarray) -> float:
    r, p, _ = pr_curve_points(y_true, scores)
    idx = np.argsort(r)
    return float(np.trapezoid(p[idx], r[idx]))


def permutation_test_auc(y_true: np.ndarray, scores: np.ndarray, n_perm: int, rng: np.random.Generator) -> Tuple[float, float, np.ndarray]:
    fpr, tpr, _ = roc_curve_points(y_true, scores)
    obs = auc_trapz(fpr, tpr)

    perms = np.zeros(n_perm, dtype=float)
    y = y_true.copy().astype(int)
    for i in range(n_perm):
        rng.shuffle(y)
        fpr_p, tpr_p, _ = roc_curve_points(y, scores)
        perms[i] = auc_trapz(fpr_p, tpr_p)

    dist = np.abs(perms - 0.5)
    pval = (np.sum(dist >= abs(obs - 0.5)) + 1) / (n_perm + 1)
    return obs, float(pval), perms


def bootstrap_ci_metric(y_true: np.ndarray, scores: np.ndarray, thr: float, n_boot: int,
                        rng: np.random.Generator) -> Dict[str, Tuple[float, float]]:
    y_true = y_true.astype(int)
    n = len(y_true)
    stats = {"auc": [], "ap": [], "f1": [], "recall": [], "specificity": [], "precision": []}
    idx = np.arange(n)

    for _ in range(n_boot):
        b = rng.choice(idx, size=n, replace=True)
        yb = y_true[b]
        sb = scores[b]
        if len(np.unique(yb[np.isfinite(sb)])) < 2:
            continue
        fpr, tpr, _ = roc_curve_points(yb, sb)
        stats["auc"].append(auc_trapz(fpr, tpr))
        stats["ap"].append(average_precision(yb, sb))
        perf = perf_from_threshold(yb, sb, thr, direction="ge")
        stats["f1"].append(perf.f1)
        stats["recall"].append(perf.recall)
        stats["specificity"].append(perf.specificity)
        stats["precision"].append(perf.precision)

    out = {}
    for k, vals in stats.items():
        if not vals:
            out[k] = (np.nan, np.nan)
        else:
            lo, hi = np.percentile(vals, [2.5, 97.5])
            out[k] = (float(lo), float(hi))
    return out


def beta_posterior(a0: float, b0: float, successes: int, failures: int) -> Tuple[float, float]:
    return a0 + successes, b0 + failures


def beta_pdf_grid(a: float, b: float, grid: np.ndarray) -> np.ndarray:
    import math
    lgamma = math.lgamma
    logB = lgamma(a) + lgamma(b) - lgamma(a + b)
    x = np.clip(grid, 1e-12, 1 - 1e-12)
    return np.exp((a - 1) * np.log(x) + (b - 1) * np.log(1 - x) - logB)


# -----------------------------
# combined models
# -----------------------------
@dataclass
class Rule:
    name: str
    metrics: List[str]
    thresholds: Dict[str, float]
    requirement: str  # "all" or "kofn"
    k: int = 0

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        rows = []
        for m in self.metrics:
            x = df[m].to_numpy(dtype=float)
            thr = float(self.thresholds[m])
            rows.append((x >= thr).astype(int))
        if not rows:
            return np.zeros(len(df), dtype=int)
        M = np.vstack(rows)
        if self.requirement == "all":
            return (M.sum(axis=0) == len(self.metrics)).astype(int)
        if self.requirement == "kofn":
            return (M.sum(axis=0) >= int(self.k)).astype(int)
        raise ValueError(f"Unknown requirement: {self.requirement}")

    def score_as_float(self, df: pd.DataFrame) -> np.ndarray:
        rows = []
        for m in self.metrics:
            x = df[m].to_numpy(dtype=float)
            thr = float(self.thresholds[m])
            rows.append((x >= thr).astype(float))
        if not rows:
            return np.zeros(len(df), dtype=float)
        return np.vstack(rows).mean(axis=0)


def pick_candidate_thresholds(x: np.ndarray, n_grid: int = 40) -> np.ndarray:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.array([np.nan])
    qs = np.linspace(0.0, 1.0, n_grid)
    thrs = np.unique(np.quantile(x, qs))
    return thrs


def choose_best_single_metric(df: pd.DataFrame, y_true: np.ndarray, metric: str) -> Tuple[Perf, float, float, float]:
    x = df[metric].to_numpy(dtype=float)
    thrs = pick_candidate_thresholds(x, n_grid=60)

    best: Optional[Perf] = None
    for thr in thrs:
        perf = perf_from_threshold(y_true, x, thr, "ge")
        if best is None:
            best = perf
            continue
        if (perf.f1 > best.f1) or (perf.f1 == best.f1 and perf.recall > best.recall) or \
           (perf.f1 == best.f1 and perf.recall == best.recall and perf.specificity > best.specificity):
            best = perf

    fpr, tpr, _ = roc_curve_points(y_true, x)
    auc = auc_trapz(fpr, tpr)
    ap = average_precision(y_true, x)
    return best, float(auc), float(ap), float(best.threshold)


def choose_best_rule(df: pd.DataFrame, y_true: np.ndarray, metrics: List[str],
                     requirement: str, k: int = 0) -> Tuple[Rule, Perf, float, float]:
    grids = {}
    for m in metrics:
        grids[m] = pick_candidate_thresholds(df[m].to_numpy(dtype=float), n_grid=30)

    keys = list(metrics)
    arrays = [grids[k] for k in keys]
    if any((a.size == 0 or (a.size == 1 and not np.isfinite(a[0]))) for a in arrays):
        r = Rule(name="invalid", metrics=metrics, thresholds={m: np.nan for m in metrics}, requirement=requirement, k=k)
        return r, Perf(np.nan, 0, 0, 0, 0, 0, 0, 0, 0), np.nan, np.nan

    best_rule: Optional[Rule] = None
    best_perf: Optional[Perf] = None

    lens = [len(a) for a in arrays]
    total = int(np.prod(lens))
    for flat in range(total):
        idxs = []
        rem = flat
        for L in reversed(lens):
            idxs.append(rem % L)
            rem //= L
        idxs = list(reversed(idxs))
        thresholds = {keys[i]: float(arrays[i][idxs[i]]) for i in range(len(keys))}
        rule = Rule(
            name=f"{requirement.upper()}_{len(metrics)}m",
            metrics=metrics,
            thresholds=thresholds,
            requirement=requirement,
            k=k,
        )
        y_pred = rule.predict(df)
        tp, tn, fp, fn = confusion_from_pred(y_true, y_pred)
        prec = safe_div(tp, tp + fp)
        rec = safe_div(tp, tp + fn)
        spec = safe_div(tn, tn + fp)
        f1 = safe_div(2 * prec * rec, prec + rec) if (prec + rec) else 0.0
        perf = Perf(threshold=np.nan, tp=tp, tn=tn, fp=fp, fn=fn, precision=prec, recall=rec, specificity=spec, f1=f1)

        if best_perf is None:
            best_perf = perf
            best_rule = rule
            continue

        if (perf.f1 > best_perf.f1) or \
           (perf.f1 == best_perf.f1 and perf.recall > best_perf.recall) or \
           (perf.f1 == best_perf.f1 and perf.recall == best_perf.recall and perf.specificity > best_perf.specificity) or \
           (perf.f1 == best_perf.f1 and perf.recall == best_perf.recall and perf.specificity == best_perf.specificity and perf.fp < best_perf.fp):
            best_perf = perf
            best_rule = rule

    score = best_rule.score_as_float(df) if best_rule else np.zeros(len(df))
    fpr, tpr, _ = roc_curve_points(y_true, score)
    auc = auc_trapz(fpr, tpr)
    ap = average_precision(y_true, score)
    return best_rule, best_perf, float(auc), float(ap)


# -----------------------------
# plotting helpers (matplotlib-only)
# -----------------------------
def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=320)
    plt.close()


def plot_metric_distributions(df: pd.DataFrame, y: np.ndarray, metric: str, out: Path) -> None:
    x = df[metric].to_numpy(dtype=float)
    pos = x[y == 1]
    neg = x[y == 0]

    plt.figure(figsize=(7, 4))
    plt.hist(neg[np.isfinite(neg)], bins=40, alpha=0.7, label="Truth negative")
    plt.hist(pos[np.isfinite(pos)], bins=40, alpha=0.7, label="Truth positive")
    plt.xlabel(metric)
    plt.ylabel("Count")
    plt.title(f"Distribution of {metric} by truth class")
    plt.legend()
    savefig(out)


def plot_roc_pr(y: np.ndarray, scores: np.ndarray, title: str, out_roc: Path, out_pr: Path) -> Tuple[float, float]:
    fpr, tpr, _ = roc_curve_points(y, scores)
    auc = auc_trapz(fpr, tpr)

    plt.figure(figsize=(5.5, 5))
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate (1 - specificity)")
    plt.ylabel("True Positive Rate (recall)")
    plt.title(f"{title}\nROC AUC = {auc:.3f}")
    savefig(out_roc)

    r, p, _ = pr_curve_points(y, scores)
    ap = average_precision(y, scores)

    plt.figure(figsize=(5.5, 5))
    plt.plot(r, p)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{title}\nAverage Precision = {ap:.3f}")
    savefig(out_pr)

    return auc, ap


def plot_confusion(tp: int, tn: int, fp: int, fn: int, title: str, out: Path) -> None:
    mat = np.array([[tp, fn],
                    [fp, tn]], dtype=int)
    plt.figure(figsize=(4.8, 4.2))
    plt.imshow(mat, interpolation="nearest")
    plt.xticks([0, 1], ["Pred +", "Pred -"])
    plt.yticks([0, 1], ["True +", "True -"])
    plt.title(title)
    for (i, j), v in np.ndenumerate(mat):
        plt.text(j, i, str(v), ha="center", va="center", fontsize=12)
    plt.colorbar(label="Count")
    plt.tight_layout()
    plt.savefig(out, dpi=320)
    plt.close()


def plot_permutation_hist(perms: np.ndarray, obs: float, pval: float, title: str, out: Path) -> None:
    plt.figure(figsize=(6.5, 4))
    plt.hist(perms, bins=40, alpha=0.75)
    plt.axvline(obs, linestyle="--")
    plt.axvline(0.5, linestyle=":")
    plt.xlabel("Permuted ROC AUC")
    plt.ylabel("Count")
    plt.title(f"{title}\nObserved AUC={obs:.3f}  Permutation p={pval:.3g}")
    savefig(out)


def plot_bootstrap_cis(ci: Dict[str, Tuple[float, float]], point: Dict[str, float], title: str, out: Path) -> None:
    keys = ["auc", "ap", "f1", "recall", "specificity", "precision"]
    labels = [k.upper() for k in keys]
    ys = np.arange(len(keys))[::-1]

    plt.figure(figsize=(7.5, 4.5))
    for yv, k in zip(ys, keys):
        lo, hi = ci.get(k, (np.nan, np.nan))
        mid = point.get(k, np.nan)
        plt.plot([lo, hi], [yv, yv], linewidth=3)
        plt.plot([mid], [yv], marker="o")
    plt.yticks(ys, labels)
    plt.xlabel("Metric value")
    plt.title(title + "\n(bootstrap 95% CI bars; point = estimate)")
    savefig(out)


def plot_beta_posteriors(post_params: Dict[str, Tuple[float, float]], out: Path, title: str) -> None:
    grid = np.linspace(0.001, 0.999, 500)
    plt.figure(figsize=(7, 4.5))
    for name, (a, b) in post_params.items():
        pdf = beta_pdf_grid(a, b, grid)
        plt.plot(grid, pdf, label=f"{name} (mean={a/(a+b):.2f})")
    plt.xlabel("Probability")
    plt.ylabel("Density")
    plt.title(title)
    plt.legend()
    savefig(out)


def plot_corr_heatmap(df: pd.DataFrame, metrics: List[str], out: Path) -> None:
    sub = df[metrics].copy()
    corr = sub.corr(numeric_only=True)
    plt.figure(figsize=(0.8 * len(metrics) + 2, 0.8 * len(metrics) + 2))
    plt.imshow(corr.values, vmin=-1, vmax=1)
    plt.colorbar(label="Pearson r")
    plt.xticks(range(len(metrics)), metrics, rotation=45, ha="right")
    plt.yticks(range(len(metrics)), metrics)
    for (i, j), v in np.ndenumerate(corr.values):
        if np.isfinite(v):
            plt.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9)
    plt.title("Correlation between candidate metrics")
    savefig(out)


# -----------------------------
# main
# -----------------------------
def build_modelling_table(truth: pd.DataFrame, bam: pd.DataFrame, genes: pd.DataFrame, blast: pd.DataFrame) -> pd.DataFrame:
    df = truth.copy()
    for other in [bam, genes, blast]:
        if other is not None and not other.empty:
            df = df.merge(other, on="base", how="left")

    if "cmv_cov_bases_mapq" in df.columns and "cmv_ref_bases" in df.columns:
        df["cmv_breadth_from_cov"] = df["cmv_cov_bases_mapq"] / df["cmv_ref_bases"].replace({0: np.nan})

    if "cmv_mapped_reads_mapq" in df.columns:
        df["log10_mapped_hq"] = np.log10(df["cmv_mapped_reads_mapq"].replace({0: np.nan}))

    if "gene_mean_depth" in df.columns:
        df["log10_gene_mean_depth"] = np.log10(df["gene_mean_depth"].replace({0: np.nan}))

    if "cmv_rpm_mapq" in df.columns:
        df["log10_rpm_hq"] = np.log10(df["cmv_rpm_mapq"].replace({0: np.nan}))

    df["y_true"] = df["y_true"].astype(int)
    return df


def run(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    ensure_dir(figdir)

    csv_dirs = [Path(p) for p in args.csv_dir]
    truth = load_truth(Path(args.truth))

    bam = load_alignment_summaries(csv_dirs)
    genes = load_gene_metrics(csv_dirs)
    blast = load_blast_top5(csv_dirs)

    df = build_modelling_table(truth, bam, genes, blast)

    candidate_metrics = [
        "cmv_mapped_reads_mapq",
        "cmv_rpm_mapq",
        "cmv_breadth_mapq",
        "cmv_breadth_from_cov",
        "gene_mean_depth",
        "gene_min_depth",
        "UL97_mean_depth",
        "UL54_mean_depth",
        "blast_cmv_mean_identity_top5",
        "blast_cmv_hits_top5",
        "log10_mapped_hq",
        "log10_rpm_hq",
        "log10_gene_mean_depth",
    ]
    candidate_metrics = [m for m in candidate_metrics if m in df.columns]

    if not candidate_metrics:
        die("No candidate metrics were found after merging. Check your csv_dir contents.")

    df = df[np.isfinite(df["y_true"].to_numpy(dtype=float))].copy()
    df.to_csv(outdir / "metrics_table.tsv", sep="\t", index=False)

    y = df["y_true"].to_numpy(dtype=int)

    plot_corr_heatmap(df, [m for m in candidate_metrics if df[m].notna().any()], figdir / "metric_correlation.png")

    rng = np.random.default_rng(args.seed)

    rows = []
    per_metric_df_rows = []
    for m in candidate_metrics:
        if not df[m].notna().any():
            continue
        best, auc, ap, thr = choose_best_single_metric(df, y, m)
        per_metric_df_rows.append({
            "metric": m,
            "best_threshold": thr,
            "tp": best.tp, "tn": best.tn, "fp": best.fp, "fn": best.fn,
            "precision": best.precision, "recall": best.recall, "specificity": best.specificity, "f1": best.f1,
            "roc_auc": auc, "avg_precision": ap,
        })

        plot_metric_distributions(df, y, m, figdir / f"dist_{m}.png")
        plot_roc_pr(y, df[m].to_numpy(dtype=float), f"Single metric: {m}", figdir / f"roc_{m}.png", figdir / f"pr_{m}.png")
        plot_confusion(best.tp, best.tn, best.fp, best.fn,
                       title=f"{m} @ thr={thr:.3g}\nF1={best.f1:.3f} Recall={best.recall:.3f} Spec={best.specificity:.3f}",
                       out=figdir / f"conf_{m}.png")

        if args.permutation_n > 0:
            obs, pval, perms = permutation_test_auc(y, df[m].to_numpy(dtype=float), args.permutation_n, rng)
            plot_permutation_hist(perms, obs, pval, f"Permutation test: {m}", figdir / f"perm_auc_{m}.png")

    per_metric_df = pd.DataFrame(per_metric_df_rows).sort_values(["f1", "recall", "specificity"], ascending=False)
    per_metric_df.to_csv(outdir / "per_metric_performance.tsv", sep="\t", index=False)

    # pool metrics for combination search
    top_f1 = per_metric_df.sort_values("f1", ascending=False).head(args.top_k_metrics)["metric"].tolist()
    top_auc = per_metric_df.sort_values("roc_auc", ascending=False).head(args.top_k_metrics)["metric"].tolist()
    pool = []
    for m in (top_f1 + top_auc):
        if m not in pool:
            pool.append(m)
    pool = pool[: max(args.top_k_metrics, 6)]

    preferred = ["cmv_mapped_reads_mapq", "cmv_breadth_mapq", "gene_mean_depth", "blast_cmv_mean_identity_top5"]
    preferred = [m for m in preferred if m in df.columns and df[m].notna().any()]
    for m in preferred:
        if m not in pool:
            pool.append(m)

    def combos(items: List[str], r: int) -> List[List[str]]:
        out = []
        n = len(items)
        idx = list(range(r))
        while True:
            out.append([items[i] for i in idx])
            for i in reversed(range(r)):
                if idx[i] != i + n - r:
                    break
            else:
                return out
            idx[i] += 1
            for j in range(i + 1, r):
                idx[j] = idx[j - 1] + 1

    model_rows = []
    best_overall = None  # (rule, perf, auc, ap)

    max_r = min(args.max_metrics, len(pool))
    for r in range(2, max_r + 1):
        for mets in combos(pool, r):
            # AND
            rule, perf, auc, ap = choose_best_rule(df, y, mets, requirement="all", k=0)
            model_rows.append({
                "model": f"AND({','.join(mets)})",
                "n_metrics": r,
                "requirement": "all",
                "k": r,
                "tp": perf.tp, "tn": perf.tn, "fp": perf.fp, "fn": perf.fn,
                "precision": perf.precision, "recall": perf.recall, "specificity": perf.specificity, "f1": perf.f1,
                "roc_auc": auc, "avg_precision": ap,
                "thresholds_json": json.dumps(rule.thresholds),
            })
            if best_overall is None:
                best_overall = (rule, perf, auc, ap)
            else:
                best_rule, best_perf, _, _ = best_overall
                if (perf.f1 > best_perf.f1) or \
                   (perf.f1 == best_perf.f1 and perf.recall > best_perf.recall) or \
                   (perf.f1 == best_perf.f1 and perf.recall == best_perf.recall and perf.specificity > best_perf.specificity) or \
                   (perf.f1 == best_perf.f1 and perf.recall == best_perf.recall and perf.specificity == best_perf.specificity and r < len(best_rule.metrics)):
                    best_overall = (rule, perf, auc, ap)

            # k-of-n (for r>=3)
            if r >= 3:
                k = max(2, r - 1)
                rule2, perf2, auc2, ap2 = choose_best_rule(df, y, mets, requirement="kofn", k=k)
                model_rows.append({
                    "model": f"{k}of{r}({','.join(mets)})",
                    "n_metrics": r,
                    "requirement": "kofn",
                    "k": k,
                    "tp": perf2.tp, "tn": perf2.tn, "fp": perf2.fp, "fn": perf2.fn,
                    "precision": perf2.precision, "recall": perf2.recall, "specificity": perf2.specificity, "f1": perf2.f1,
                    "roc_auc": auc2, "avg_precision": ap2,
                    "thresholds_json": json.dumps(rule2.thresholds),
                })
                best_rule, best_perf, _, _ = best_overall
                if (perf2.f1 > best_perf.f1) or \
                   (perf2.f1 == best_perf.f1 and perf2.recall > best_perf.recall) or \
                   (perf2.f1 == best_perf.f1 and perf2.recall == best_perf.recall and perf2.specificity > best_perf.specificity) or \
                   (perf2.f1 == best_perf.f1 and perf2.recall == best_perf.recall and perf2.specificity == best_perf.specificity and r < len(best_rule.metrics)):
                    best_overall = (rule2, perf2, auc2, ap2)

    model_df = pd.DataFrame(model_rows).sort_values(["f1", "recall", "specificity", "n_metrics"],
                                                    ascending=[False, False, False, True])
    model_df.to_csv(outdir / "model_comparison.tsv", sep="\t", index=False)

    if best_overall is None:
        die("No combined model could be fit (check that metrics have data).")

    best_rule, best_perf, best_auc, best_ap = best_overall
        # -----------------------------
    # EXTRA TABLES: show how final exact thresholds were chosen
    # -----------------------------
    def evaluate_rule_with_thresholds(thr_map: Dict[str, float]) -> Dict[str, float]:
        """Evaluate the FINAL rule structure (same metrics, same requirement/k) using a given threshold dict."""
        tmp_rule = Rule(
            name="FINAL_RULE_REEVAL",
            metrics=best_rule.metrics,
            thresholds=thr_map,
            requirement=best_rule.requirement,
            k=best_rule.k,
        )
        yp = tmp_rule.predict(df)
        tp_, tn_, fp_, fn_ = confusion_from_pred(y, yp)
        prec_ = safe_div(tp_, tp_ + fp_)
        rec_ = safe_div(tp_, tp_ + fn_)
        spec_ = safe_div(tn_, tn_ + fp_)
        f1_ = safe_div(2 * prec_ * rec_, prec_ + rec_) if (prec_ + rec_) else 0.0
        score_ = tmp_rule.score_as_float(df)
        fpr_, tpr_, _ = roc_curve_points(y, score_)
        auc_ = auc_trapz(fpr_, tpr_)
        ap_ = average_precision(y, score_)
        return {
            "tp": tp_, "tn": tn_, "fp": fp_, "fn": fn_,
            "precision": prec_, "recall": rec_, "specificity": spec_, "f1": f1_,
            "roc_auc": auc_, "avg_precision": ap_,
        }

    # Recreate the EXACT grids used by choose_best_rule for the final metrics
    final_grids = {}
    for m in best_rule.metrics:
        final_grids[m] = pick_candidate_thresholds(df[m].to_numpy(dtype=float), n_grid=30)

    # Table 1: exact final thresholds + grid info
    rows_final = []
    for m in best_rule.metrics:
        grid = final_grids[m]
        final_thr = float(best_rule.thresholds[m])

        # Find closest grid value and index (should usually be exact)
        if grid.size == 0 or (grid.size == 1 and not np.isfinite(grid[0])):
            rows_final.append({
                "metric": m,
                "final_threshold": final_thr,
                "grid_size": int(grid.size),
                "closest_grid_value": np.nan,
                "closest_grid_index": np.nan,
                "difference_to_grid": np.nan,
                "grid_min": np.nan,
                "grid_max": np.nan,
            })
            continue

        idx = int(np.nanargmin(np.abs(grid - final_thr)))
        closest = float(grid[idx])
        rows_final.append({
            "metric": m,
            "final_threshold": final_thr,
            "grid_size": int(grid.size),
            "closest_grid_value": closest,
            "closest_grid_index": idx,
            "difference_to_grid": float(final_thr - closest),
            "grid_min": float(np.nanmin(grid)),
            "grid_max": float(np.nanmax(grid)),
        })

    final_thresholds_table = pd.DataFrame(rows_final)
    final_thresholds_table.to_csv(outdir / "final_thresholds_table.tsv", sep="\t", index=False)

    # Table 2: local sensitivity around each final threshold (neighbours in the grid)
    base_eval = evaluate_rule_with_thresholds(best_rule.thresholds)
    rows_local = []

    for m in best_rule.metrics:
        grid = final_grids[m]
        if grid.size < 2 or not np.isfinite(best_rule.thresholds[m]):
            continue

        final_thr = float(best_rule.thresholds[m])
        idx = int(np.nanargmin(np.abs(grid - final_thr)))

        # take neighbouring grid points around the selected one
        neighbour_idxs = sorted(set([idx - 2, idx - 1, idx, idx + 1, idx + 2]))
        neighbour_idxs = [i for i in neighbour_idxs if 0 <= i < len(grid)]

        for i in neighbour_idxs:
            thr_map = dict(best_rule.thresholds)
            thr_map[m] = float(grid[i])
            ev = evaluate_rule_with_thresholds(thr_map)

            rows_local.append({
                "metric_varied": m,
                "grid_index": i,
                "threshold_tested": float(grid[i]),
                "is_final_choice": (i == idx),
                # performance
                "tp": ev["tp"], "tn": ev["tn"], "fp": ev["fp"], "fn": ev["fn"],
                "precision": ev["precision"], "recall": ev["recall"], "specificity": ev["specificity"], "f1": ev["f1"],
                "roc_auc": ev["roc_auc"], "avg_precision": ev["avg_precision"],
                # deltas vs final
                "delta_f1_vs_final": float(ev["f1"] - base_eval["f1"]),
                "delta_recall_vs_final": float(ev["recall"] - base_eval["recall"]),
                "delta_spec_vs_final": float(ev["specificity"] - base_eval["specificity"]),
                "delta_fp_vs_final": int(ev["fp"] - base_eval["fp"]),
                "delta_fn_vs_final": int(ev["fn"] - base_eval["fn"]),
            })

    local_sensitivity_table = pd.DataFrame(rows_local).sort_values(
        ["metric_varied", "grid_index"]
    )
    local_sensitivity_table.to_csv(outdir / "final_threshold_local_sensitivity.tsv", sep="\t", index=False)

    y_pred = best_rule.predict(df)
    tp, tn, fp, fn = confusion_from_pred(y, y_pred)

    plot_confusion(tp, tn, fp, fn,
                   title=f"Final combined rule: {best_rule.requirement} ({len(best_rule.metrics)} metrics)\n"
                         f"F1={best_perf.f1:.3f} Recall={best_perf.recall:.3f} Spec={best_perf.specificity:.3f}",
                   out=figdir / "final_confusion.png")

    score = best_rule.score_as_float(df)
    plot_roc_pr(y, score, "Final combined rule (score = fraction of metrics passing)",
                figdir / "final_roc.png", figdir / "final_pr.png")

    if args.permutation_n > 0:
        obs, pval, perms = permutation_test_auc(y, score, args.permutation_n, rng)
        plot_permutation_hist(perms, obs, pval, "Permutation test: final combined rule", figdir / "final_perm_auc.png")

    # Bootstrap on rule score threshold 0.5 (score >=0.5 means >= half the metrics pass)
    ci = bootstrap_ci_metric(y, score, thr=0.5, n_boot=args.bootstrap_n, rng=rng)
    point = {
        "auc": best_auc,
        "ap": best_ap,
        "f1": best_perf.f1,
        "recall": best_perf.recall,
        "specificity": best_perf.specificity,
        "precision": best_perf.precision,
    }
    plot_bootstrap_cis(ci, point, "Final model performance", figdir / "final_bootstrap_ci.png")

    a_sens, b_sens = beta_posterior(1, 1, tp, fn)
    a_spec, b_spec = beta_posterior(1, 1, tn, fp)
    a_ppv, b_ppv = beta_posterior(1, 1, tp, fp)
    a_npv, b_npv = beta_posterior(1, 1, tn, fn)
    post = {
        "Sensitivity": (a_sens, b_sens),
        "Specificity": (a_spec, b_spec),
        "PPV": (a_ppv, b_ppv),
        "NPV": (a_npv, b_npv),
    }
    plot_beta_posteriors(post, figdir / "final_bayesian_posteriors.png",
                         title="Bayesian posteriors at final threshold (Beta(1,1) prior)")

    out_json = {
        "rule_name": best_rule.name,
        "requirement": best_rule.requirement,
        "k": best_rule.k,
        "metrics": best_rule.metrics,
        "thresholds": best_rule.thresholds,
        "operating_point": {
            "tp": tp, "tn": tn, "fp": fp, "fn": fn,
            "precision": best_perf.precision,
            "recall": best_perf.recall,
            "specificity": best_perf.specificity,
            "f1": best_perf.f1,
            "roc_auc": best_auc,
            "avg_precision": best_ap,
        },
        "bootstrap_ci": {k: {"lo": v[0], "hi": v[1]} for k, v in ci.items()},
        "bayesian_posteriors": {k: {"a": float(v[0]), "b": float(v[1]), "mean": float(v[0] / (v[0] + v[1]))} for k, v in post.items()},
        "inputs": {"csv_dirs": [str(p) for p in csv_dirs], "truth": str(args.truth)},
    }
    with open(outdir / "thresholds.json", "w") as fh:
        json.dump(out_json, fh, indent=2)

    summary = pd.DataFrame([{
        "model": f"FINAL_{best_rule.requirement}_{best_rule.k}of{len(best_rule.metrics)}",
        "metrics": ",".join(best_rule.metrics),
        "thresholds": json.dumps(best_rule.thresholds),
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "precision": best_perf.precision,
        "recall": best_perf.recall,
        "specificity": best_perf.specificity,
        "f1": best_perf.f1,
        "roc_auc": best_auc,
        "avg_precision": best_ap,
    }])
    summary.to_csv(outdir / "final_model_summary.tsv", sep="\t", index=False)

    print("\n🎉 Threshold optimisation complete.")
    print(f"Outputs: {outdir}")
    print(f"Figures: {figdir}")
    print(f"Final rule JSON: {outdir / 'thresholds.json'}")


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Optimise a multi-metric CMV detection threshold + generate dissertation-ready plots."
    )
    ap.add_argument("--truth", required=True, help="Path to truth/metadata.tsv from the simulator sample set.")
    ap.add_argument("--csv-dir", required=True, nargs="+",
                    help="One or more pipeline CSV output directories (pass multiple if you ran in batches).")
    ap.add_argument("--outdir", required=True, help="Output directory for tables + figures.")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-metrics", type=int, default=4,
                    help="Max number of metrics to consider in combined rules (default 4).")
    ap.add_argument("--top-k-metrics", type=int, default=6,
                    help="How many top single metrics to pool for combinations (default 6).")
    ap.add_argument("--permutation-n", type=int, default=2000,
                    help="Number of permutations for AUC permutation tests (0 disables).")
    ap.add_argument("--bootstrap-n", type=int, default=2000,
                    help="Number of bootstrap resamples for CI plots.")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    run(args)


if __name__ == "__main__":
    main()