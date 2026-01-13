#!/usr/bin/env python3
"""
CMV alignment-threshold optimisation (1D + 2D) with:
  - ROC + AUC (rank/Mann–Whitney, no sklearn)
  - PR curve + Average Precision (AP, no sklearn)
  - Permutation test for AUC (optional)
  - Bootstrap CIs for AUC, AP, selected threshold, and sens/spec/prec/F1 at that threshold
  - 2D threshold search: mapped_hq + (breadth OR covered_bases) with heatmap + best region output

Notes:
  - If your CMV BAM contains ONLY mapped reads, RPM computed as mapped/total_in_bam can be misleading.
    The script prints a warning if it detects many samples with total_reads_in_cmv_bam ~ mapped_reads.
"""

import argparse
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------
# sample name canonicalization (matches your style)
# ----------------------------
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


# ----------------------------
# basic metrics helpers
# ----------------------------
def confusion_counts(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn


def safe_div(a, b):
    return float(a) / float(b) if b else float("nan")


def metrics_from_counts(tp, tn, fp, fn):
    sens = safe_div(tp, tp + fn)          # recall / TPR
    spec = safe_div(tn, tn + fp)          # TNR
    prec = safe_div(tp, tp + fp)          # PPV
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    fpr = safe_div(fp, fp + tn)
    acc = safe_div(tp + tn, tp + tn + fp + fn)
    youden = sens + spec - 1 if (not np.isnan(sens) and not np.isnan(spec)) else float("nan")
    return sens, spec, prec, f1, fpr, acc, youden


# ----------------------------
# ROC / AUC (no sklearn)
# ----------------------------
def roc_auc_rank(y_true: np.ndarray, scores: np.ndarray) -> float:
    """AUC via rank/Mann–Whitney. Robust, no sklearn."""
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)

    pos = scores[y_true == 1]
    neg = scores[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")

    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)

    # average ranks for ties
    uniq, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    for uidx, cnt in enumerate(counts):
        if cnt > 1:
            idxs = np.where(inv == uidx)[0]
            ranks[idxs] = ranks[idxs].mean()

    sum_ranks_pos = ranks[y_true == 1].sum()
    n_pos = len(pos)
    n_neg = len(neg)
    auc = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return float(auc)


def roc_points_from_scores(y_true: np.ndarray, scores: np.ndarray):
    """ROC curve by thresholding scores (>= t). Returns array with columns: FPR, TPR, threshold."""
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)

    ts = np.unique(scores)
    ts = np.sort(ts)
    ts = np.concatenate(([ts[0] - 1e-9], ts, [ts[-1] + 1e-9]))

    pts = []
    for t in ts:
        y_pred = (scores >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y_true, y_pred)
        sens, spec, prec, f1, fpr, acc, youden = metrics_from_counts(tp, tn, fp, fn)
        pts.append((fpr, sens, t))
    pts = np.array(pts, dtype=float)
    pts = pts[np.argsort(pts[:, 0])]
    return pts


# ----------------------------
# PR curve + Average Precision (no sklearn)
# ----------------------------
def pr_points_from_scores(y_true: np.ndarray, scores: np.ndarray):
    """
    PR curve by thresholding scores (>= t).
    Returns arrays: recall, precision, threshold sorted by recall ascending for plotting.
    """
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)

    ts = np.unique(scores)
    ts = np.sort(ts)
    ts = np.concatenate(([ts[0] - 1e-9], ts, [ts[-1] + 1e-9]))

    rows = []
    for t in ts:
        y_pred = (scores >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y_true, y_pred)
        sens, spec, prec, f1, fpr, acc, youden = metrics_from_counts(tp, tn, fp, fn)
        # sens = recall
        rows.append((sens, prec, t))
    arr = np.array(rows, dtype=float)

    # For PR plotting it's common to sort by recall
    arr = arr[np.argsort(arr[:, 0])]
    recall = arr[:, 0]
    precision = arr[:, 1]
    threshold = arr[:, 2]
    return recall, precision, threshold


def average_precision(y_true: np.ndarray, scores: np.ndarray) -> float:
    """
    Average Precision (AP) computed like sklearn's definition:
      sort by score desc, then sum over recall increases: AP = Σ (Δrecall * precision_at_step)
    """
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)

    n_pos = int((y_true == 1).sum())
    if n_pos == 0:
        return float("nan")

    order = np.argsort(-scores)  # descending
    y_sorted = y_true[order]

    tp = 0
    fp = 0
    prev_recall = 0.0
    ap = 0.0

    # Walk down ranked list; every item is a threshold step
    for i in range(len(y_sorted)):
        if y_sorted[i] == 1:
            tp += 1
        else:
            fp += 1

        recall = tp / n_pos
        precision = tp / (tp + fp) if (tp + fp) else 1.0

        # Only add when recall increases (i.e. at positives)
        if y_sorted[i] == 1:
            ap += (recall - prev_recall) * precision
            prev_recall = recall

    return float(ap)


# ----------------------------
# samtools helpers
# ----------------------------
def run_cmd(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n  {' '.join(cmd)}\n\nSTDERR:\n{p.stderr}")
    return p.stdout


def samtools_count(bam: Path, mapped_only=False, mapq=None) -> int:
    cmd = ["samtools", "view", "-c"]
    if mapped_only:
        cmd += ["-F", "4"]
    if mapq is not None:
        cmd += ["-q", str(int(mapq))]
    cmd += [str(bam)]
    out = run_cmd(cmd).strip()
    return int(out) if out else 0


def samtools_idxstats_lengths_and_mapped(bam: Path):
    out = run_cmd(["samtools", "idxstats", str(bam)])
    genome_len = 0
    mapped_sum = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        contig, length, mapped = parts[0], parts[1], parts[2]
        if contig == "*":
            continue
        try:
            L = int(length)
            M = int(mapped)
        except ValueError:
            continue
        genome_len += L
        mapped_sum += M
    return genome_len, mapped_sum


def breadth_of_coverage(bam: Path, mapq: int, depth_min: int = 1, method: str = "coverage"):
    """
    Breadth = covered_bases / reference_length.
    method:
      - "coverage": uses `samtools coverage` if available (fast)
      - "depth": uses `samtools depth -a` (slower but widely available)
    """
    genome_len, _ = samtools_idxstats_lengths_and_mapped(bam)
    if genome_len <= 0:
        return float("nan"), 0, 0

    if method == "coverage":
        try:
            out = run_cmd(["samtools", "coverage", "-q", str(int(mapq)), str(bam)])
            covbases_sum = 0
            len_sum = 0
            for line in out.splitlines():
                if not line or line.startswith("#") or line.lower().startswith("rname"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                rname = parts[0]
                if rname == "*":
                    continue
                start = int(parts[1])
                end = int(parts[2])
                covbases = int(parts[4])
                length = end - start + 1
                covbases_sum += covbases
                len_sum += length
            if len_sum > 0:
                return covbases_sum / len_sum, covbases_sum, len_sum
        except Exception:
            method = "depth"

    out = run_cmd(["samtools", "depth", "-a", "-q", str(int(mapq)), str(bam)])
    covered = 0
    total = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            d = int(parts[2])
        except ValueError:
            continue
        total += 1
        if d >= depth_min:
            covered += 1

    if total == 0:
        return float("nan"), 0, 0
    return covered / total, covered, total


# ----------------------------
# IO helpers
# ----------------------------
def load_truth(truth_path: Path) -> pd.DataFrame:
    """
    Needs columns:
      sample, truth_cmv_present  (or cmv_present)
    """
    df = pd.read_csv(truth_path, sep="\t")
    if "sample" not in df.columns:
        raise ValueError("Truth TSV must contain a 'sample' column.")
    df["base"] = df["sample"].astype(str).apply(canonical_sample)

    if "truth_cmv_present" in df.columns:
        y = df["truth_cmv_present"]
    elif "cmv_present" in df.columns:
        y = df["cmv_present"]
    else:
        raise ValueError("Truth TSV must contain 'truth_cmv_present' or 'cmv_present' column.")

    df["truth_cmv_present"] = pd.to_numeric(y, errors="coerce").fillna(0).astype(int)
    df = df.sort_values("base").drop_duplicates("base", keep="first").reset_index(drop=True)
    return df[["base", "truth_cmv_present"]]


def find_bams(cmv_bam_dir: Path):
    bams = sorted([p for p in cmv_bam_dir.glob("*.bam") if not str(p).endswith(".bai")])
    if not bams:
        raise FileNotFoundError(f"No .bam files found in: {cmv_bam_dir}")
    return bams


def pick_score(df: pd.DataFrame, metric: str):
    if metric == "mapped_hq":
        return df["cmv_mapped_reads_mapq"].to_numpy(dtype=float)
    if metric == "rpm_hq":
        return df["cmv_rpm_mapq"].to_numpy(dtype=float)
    if metric == "breadth":
        return df["cmv_breadth_mapq"].to_numpy(dtype=float)
    if metric == "cov_bases":
        return df["cmv_cov_bases_mapq"].to_numpy(dtype=float)
    if metric == "combo":
        # simple combined score: log1p(mapped) + 10*breadth
        return np.log1p(df["cmv_mapped_reads_mapq"].to_numpy(dtype=float)) + 10.0 * df["cmv_breadth_mapq"].to_numpy(dtype=float)
    raise ValueError(f"Unknown metric: {metric}")


# ----------------------------
# threshold selection (1D)
# ----------------------------
def threshold_scan_1d(y: np.ndarray, scores: np.ndarray):
    thresholds = np.unique(scores)
    thresholds = np.sort(np.concatenate([thresholds, thresholds + 1e-12]))
    rows = []
    for t in thresholds:
        y_pred = (scores >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y, y_pred)
        sens, spec, prec, f1, fpr, acc, youden = metrics_from_counts(tp, tn, fp, fn)
        rows.append({
            "threshold_score": float(t),
            "TP": tp, "TN": tn, "FP": fp, "FN": fn,
            "sensitivity": sens,
            "specificity": spec,
            "precision": prec,
            "F1": f1,
            "FPR": fpr,
            "accuracy": acc,
            "youden_J": youden,
        })
    perf = pd.DataFrame(rows).sort_values("threshold_score").reset_index(drop=True)
    return perf


def select_threshold_1d(perf: pd.DataFrame, rule: str, min_sensitivity: float, max_fpr: float):
    best = None
    if rule == "youden":
        # Maximise Youden; tie-breaker: prefer LOWER threshold (more permissive) if identical
        perf2 = perf.copy()
        perf2["youden_J_filled"] = perf2["youden_J"].fillna(-1e18)
        best_score = perf2["youden_J_filled"].max()
        candidates = perf2[perf2["youden_J_filled"] == best_score]
        if not candidates.empty:
            best = candidates.sort_values("threshold_score", ascending=True).iloc[0].to_dict()
        return best

    if rule == "minscore_sens_fpr":
        feas = perf[(perf["sensitivity"] >= min_sensitivity) & (perf["FPR"] <= max_fpr)]
        if feas.empty:
            return None
        # choose smallest threshold meeting constraints
        best = feas.sort_values("threshold_score", ascending=True).iloc[0].to_dict()
        return best

    raise ValueError(f"Unknown rule: {rule}")


# ----------------------------
# bootstrap CIs
# ----------------------------
def percentile_ci(x, alpha=0.05):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (float("nan"), float("nan"))
    lo = np.percentile(x, 100 * (alpha / 2))
    hi = np.percentile(x, 100 * (1 - alpha / 2))
    return float(lo), float(hi)


def bootstrap_summary(y, scores, rule, min_sens, max_fpr, n_boot=2000, seed=1):
    """
    Returns dict with bootstrap arrays and CIs for:
      - AUC
      - AP
      - selected threshold
      - sens/spec/prec/F1 at selected threshold
    """
    rng = np.random.default_rng(seed)
    n = len(y)

    aucs, aps, ths = [], [], []
    senss, specs, precs, f1s, fprs = [], [], [], [], []
    skipped = 0

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        sb = scores[idx]

        # Need at least one pos and one neg for ROC AUC
        if (yb == 1).sum() == 0 or (yb == 0).sum() == 0:
            skipped += 1
            continue

        aucs.append(roc_auc_rank(yb, sb))
        aps.append(average_precision(yb, sb))

        perf = threshold_scan_1d(yb, sb)
        best = select_threshold_1d(perf, rule, min_sens, max_fpr)
        if best is None:
            ths.append(float("nan"))
            senss.append(float("nan"))
            specs.append(float("nan"))
            precs.append(float("nan"))
            f1s.append(float("nan"))
            fprs.append(float("nan"))
        else:
            ths.append(float(best["threshold_score"]))
            senss.append(float(best["sensitivity"]))
            specs.append(float(best["specificity"]))
            precs.append(float(best["precision"]))
            f1s.append(float(best["F1"]))
            fprs.append(float(best["FPR"]))

    out = {
        "aucs": np.asarray(aucs, dtype=float),
        "aps": np.asarray(aps, dtype=float),
        "ths": np.asarray(ths, dtype=float),
        "sens": np.asarray(senss, dtype=float),
        "spec": np.asarray(specs, dtype=float),
        "prec": np.asarray(precs, dtype=float),
        "f1": np.asarray(f1s, dtype=float),
        "fpr": np.asarray(fprs, dtype=float),
        "skipped": skipped,
    }

    out["auc_ci95"] = percentile_ci(out["aucs"])
    out["ap_ci95"] = percentile_ci(out["aps"])
    out["th_ci95"] = percentile_ci(out["ths"])
    out["sens_ci95"] = percentile_ci(out["sens"])
    out["spec_ci95"] = percentile_ci(out["spec"])
    out["prec_ci95"] = percentile_ci(out["prec"])
    out["f1_ci95"] = percentile_ci(out["f1"])
    out["fpr_ci95"] = percentile_ci(out["fpr"])
    return out


# ----------------------------
# 2D threshold search
# ----------------------------
def threshold_scan_2d(y, mapped_hq, second, second_name: str):
    """
    Scan all threshold pairs:
      pred = (mapped_hq >= Tm) AND (second >= Ts)
    Returns dataframe with Tm, Ts and metrics.
    """
    y = np.asarray(y).astype(int)
    mapped_hq = np.asarray(mapped_hq).astype(float)
    second = np.asarray(second).astype(float)

    t_m = np.sort(np.unique(mapped_hq))
    t_s = np.sort(np.unique(second))

    rows = []
    for tm in t_m:
        for ts in t_s:
            y_pred = ((mapped_hq >= tm) & (second >= ts)).astype(int)
            tp, tn, fp, fn = confusion_counts(y, y_pred)
            sens, spec, prec, f1, fpr, acc, youden = metrics_from_counts(tp, tn, fp, fn)
            rows.append({
                "Tm_mapped_hq": float(tm),
                f"Ts_{second_name}": float(ts),
                "TP": tp, "TN": tn, "FP": fp, "FN": fn,
                "sensitivity": sens,
                "specificity": spec,
                "precision": prec,
                "F1": f1,
                "FPR": fpr,
                "accuracy": acc,
                "youden_J": youden,
            })
    df = pd.DataFrame(rows)
    return df, t_m, t_s


def select_threshold_2d(scan2d: pd.DataFrame, rule: str, min_sens: float, max_fpr: float, second_name: str):
    if scan2d.empty:
        return None

    if rule == "youden":
        s = scan2d["youden_J"].fillna(-1e18)
        best = scan2d.loc[s.idxmax()].to_dict()
        return best

    if rule == "max_f1":
        s = scan2d["F1"].fillna(-1e18)
        best = scan2d.loc[s.idxmax()].to_dict()
        return best

    if rule == "minscore_sens_fpr":
        feas = scan2d[(scan2d["sensitivity"] >= min_sens) & (scan2d["FPR"] <= max_fpr)]
        if feas.empty:
            return None
        # lexicographic: smallest mapped threshold, then smallest second threshold
        cols = ["Tm_mapped_hq", f"Ts_{second_name}"]
        best = feas.sort_values(cols, ascending=True).iloc[0].to_dict()
        return best

    raise ValueError(f"Unknown 2D rule: {rule}")


# ----------------------------
# plotting helpers (well-labeled)
# ----------------------------
def save_roc_plot(path: Path, roc_pts, auc, auc_ci, metric, n_pos, n_neg, p_perm=None):
    plt.figure(figsize=(6.2, 5.2))
    plt.plot(roc_pts[:, 0], roc_pts[:, 1], label="ROC curve")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Chance (AUC=0.5)")
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    title = f"ROC curve (score = {metric})\nAUC = {auc:.3f} (95% CI {auc_ci[0]:.3f}–{auc_ci[1]:.3f})"
    if p_perm is not None and np.isfinite(p_perm):
        title += f" | Permutation p = {p_perm:.3g}"
    plt.title(title)
    plt.legend(loc="lower right", frameon=True)
    plt.text(0.02, 0.02, f"n_pos={n_pos}, n_neg={n_neg}", transform=plt.gca().transAxes)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def save_pr_plot(path: Path, recall, precision, ap, ap_ci, metric, prevalence):
    plt.figure(figsize=(6.2, 5.2))
    plt.plot(recall, precision, label="PR curve")
    # baseline precision = prevalence
    plt.axhline(prevalence, linestyle="--", label=f"Baseline (prevalence={prevalence:.3f})")
    plt.xlabel("Recall (Sensitivity)")
    plt.ylabel("Precision (PPV)")
    plt.ylim(0, 1.02)
    plt.xlim(0, 1.0)
    plt.title(f"Precision–Recall curve (score = {metric})\nAP = {ap:.3f} (95% CI {ap_ci[0]:.3f}–{ap_ci[1]:.3f})")
    plt.legend(loc="lower left", frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def save_metrics_vs_threshold_plot(path: Path, perf: pd.DataFrame, metric: str, best_threshold=None, rule_desc=""):
    plt.figure(figsize=(7.0, 5.2))
    plt.plot(perf["threshold_score"], perf["sensitivity"], label="Sensitivity (Recall)")
    plt.plot(perf["threshold_score"], perf["specificity"], label="Specificity")
    plt.plot(perf["threshold_score"], perf["precision"], label="Precision (PPV)")
    plt.plot(perf["threshold_score"], perf["F1"], label="F1 score")
    plt.plot(perf["threshold_score"], perf["youden_J"], label="Youden's J")

    if best_threshold is not None and np.isfinite(best_threshold):
        plt.axvline(best_threshold, linestyle="--", label=f"Selected threshold = {best_threshold:.4g}")

    plt.xlabel("Score threshold T (call positive if score ≥ T)")
    plt.ylabel("Metric value")
    ttl = f"Threshold scan (score = {metric})"
    if rule_desc:
        ttl += f"\nSelection rule: {rule_desc}"
    plt.title(ttl)
    plt.legend(loc="best", frameon=True, ncol=2)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def save_bootstrap_hist(path: Path, values: np.ndarray, xlabel: str, title: str, vline=None):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    plt.figure(figsize=(6.2, 5.2))
    if values.size:
        plt.hist(values, bins=45)
    if vline is not None and np.isfinite(vline):
        plt.axvline(vline, linestyle="--", label="Observed")
        plt.legend(frameon=True)
    plt.xlabel(xlabel)
    plt.ylabel("Bootstrap count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def save_2d_heatmap(path: Path, scan2d: pd.DataFrame, t_m: np.ndarray, t_s: np.ndarray,
                    second_name: str, value_col: str, best_pair=None, title_extra=""):
    # build matrix [len(t_m) x len(t_s)]
    mat = np.full((len(t_m), len(t_s)), np.nan, dtype=float)
    # indexing maps
    idx_m = {v: i for i, v in enumerate(t_m)}
    idx_s = {v: j for j, v in enumerate(t_s)}

    for _, row in scan2d.iterrows():
        i = idx_m[row["Tm_mapped_hq"]]
        j = idx_s[row[f"Ts_{second_name}"]]
        mat[i, j] = row[value_col]

    plt.figure(figsize=(7.4, 5.6))
    im = plt.imshow(mat, origin="lower", aspect="auto")
    plt.colorbar(im, label=value_col)

    plt.xlabel(f"Threshold Ts ({second_name})")
    plt.ylabel("Threshold Tm (mapped_hq)")
    plt.title(f"2D threshold scan heatmap: {value_col}\n"
              f"Positive if (mapped_hq ≥ Tm) AND ({second_name} ≥ Ts){title_extra}")

    # Use a sparse tick set for readability
    def sparse_ticks(vals, max_ticks=8):
        if len(vals) <= max_ticks:
            idx = np.arange(len(vals))
        else:
            idx = np.unique(np.round(np.linspace(0, len(vals) - 1, max_ticks)).astype(int))
        return idx, [f"{vals[i]:.4g}" for i in idx]

    xt_idx, xt_lbl = sparse_ticks(t_s, max_ticks=8)
    yt_idx, yt_lbl = sparse_ticks(t_m, max_ticks=8)
    plt.xticks(xt_idx, xt_lbl, rotation=45, ha="right")
    plt.yticks(yt_idx, yt_lbl)

    if best_pair is not None:
        tm = float(best_pair["Tm_mapped_hq"])
        ts = float(best_pair[f"Ts_{second_name}"])
        plt.scatter([idx_s[ts]], [idx_m[tm]], marker="x")
        plt.text(idx_s[ts], idx_m[tm], "  best", va="center")

    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


# ----------------------------
# main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Alignment-based CMV detection: ROC + PR + bootstrap CIs + optional 2D threshold search."
    )
    ap.add_argument("--truth", required=True, help="Truth TSV with sample + truth_cmv_present (0/1).")
    ap.add_argument("--cmv-bam-dir", required=True, help="Directory containing CMV-alignment BAMs.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--mapq", type=int, default=30, help="MAPQ cutoff for CMV mapped reads and breadth.")
    ap.add_argument("--breadth-depth-min", type=int, default=1, help="Depth threshold for breadth (>= this counts as covered).")
    ap.add_argument("--breadth-method", choices=["coverage", "depth"], default="coverage", help="Breadth method.")

    ap.add_argument("--metric", choices=["mapped_hq", "rpm_hq", "breadth", "cov_bases", "combo"], default="mapped_hq",
                    help="Which alignment metric to use as 1D ROC/PR score.")
    ap.add_argument("--rule", choices=["youden", "minscore_sens_fpr"], default="youden",
                    help="Threshold selection rule for 1D score.")
    ap.add_argument("--min-sensitivity", type=float, default=0.99, help="Used when rule=minscore_sens_fpr.")
    ap.add_argument("--max-fpr", type=float, default=0.05, help="Used when rule=minscore_sens_fpr.")

    ap.add_argument("--permute", type=int, default=0, help="Permutation iterations for AUC p-value (0 disables).")
    ap.add_argument("--bootstrap", type=int, default=2000, help="Bootstrap resamples for CI estimation (0 disables).")
    ap.add_argument("--seed", type=int, default=1)

    # 2D search options
    ap.add_argument("--two-dim", choices=["none", "mapped_breadth", "mapped_covbases"], default="none",
                    help="Run 2D threshold search: mapped_hq + breadth OR mapped_hq + covered_bases.")
    ap.add_argument("--rule2d", choices=["youden", "max_f1", "minscore_sens_fpr"], default="youden",
                    help="Selection rule for 2D thresholds.")
    ap.add_argument("--top2d", type=int, default=50, help="Write top-N 2D threshold pairs by Youden (or chosen objective).")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load truth
    truth = load_truth(Path(args.truth))

    # Find BAMs and compute metrics
    bams = find_bams(Path(args.cmv_bam_dir))

    rows = []
    for bam in bams:
        sample = re.sub(r"\.bam$", "", bam.name, flags=re.IGNORECASE)
        base = canonical_sample(sample)

        total_reads_bam = samtools_count(bam, mapped_only=False, mapq=None)
        mapped_hq = samtools_count(bam, mapped_only=True, mapq=args.mapq)
        rpm_hq = (mapped_hq / total_reads_bam * 1e6) if total_reads_bam else 0.0

        breadth, cov_bases, ref_bases = breadth_of_coverage(
            bam, mapq=args.mapq, depth_min=args.breadth_depth_min, method=args.breadth_method
        )

        rows.append({
            "base": base,
            "bam": str(bam),
            "total_reads_in_cmv_bam": int(total_reads_bam),
            "cmv_mapped_reads_mapq": int(mapped_hq),
            "cmv_rpm_mapq": float(rpm_hq),
            "cmv_breadth_mapq": float(breadth),
            "cmv_cov_bases_mapq": int(cov_bases),
            "cmv_ref_bases": int(ref_bases),
        })

    align_df = pd.DataFrame(rows).sort_values("base").reset_index(drop=True)

    # Basic RPM sanity warning
    if "total_reads_in_cmv_bam" in align_df.columns and "cmv_mapped_reads_mapq" in align_df.columns:
        # If BAM is mapped-only, total ~= mapped for many samples.
        # We'll warn if >30% samples have total <= mapped+2 and mapped > 0.
        near_equal = ((align_df["total_reads_in_cmv_bam"] <= (align_df["cmv_mapped_reads_mapq"] + 2)) &
                      (align_df["cmv_mapped_reads_mapq"] > 0)).mean()
        if near_equal > 0.30:
            print("⚠️ Warning: Many samples have total_reads_in_cmv_bam ~= mapped_hq.")
            print("   If your CMV BAMs are mapped-only, RPM computed using BAM total is not meaningful.")
            print("   Prefer mapped_hq + breadth/cov_bases, or provide an external denominator.")

    # Merge truth
    df = truth.merge(align_df, on="base", how="left")
    if df["bam"].isna().any():
        missing = df[df["bam"].isna()]["base"].tolist()
        print("⚠️ Warning: missing BAMs for these bases in truth table:")
        print("  " + ", ".join(missing))

    # Drop rows without BAMs for analysis
    df_roc = df.dropna(subset=["bam"]).copy()
    y = df_roc["truth_cmv_present"].to_numpy(dtype=int)

    scores = pick_score(df_roc, args.metric)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    prevalence = safe_div(n_pos, len(y))

    # 1D ROC / PR
    auc = roc_auc_rank(y, scores)
    roc_pts = roc_points_from_scores(y, scores)
    recall, precision, _thr_pr = pr_points_from_scores(y, scores)
    ap_score = average_precision(y, scores)

    # Threshold scan & selection (1D)
    perf = threshold_scan_1d(y, scores)
    best = select_threshold_1d(perf, args.rule, args.min_sensitivity, args.max_fpr)

    # Permutation test for AUC (optional)
    p_perm = None
    auc_null = None
    if args.permute and args.permute > 0:
        rng = np.random.default_rng(args.seed)
        aucs = []
        for _ in range(args.permute):
            y_perm = rng.permutation(y)
            aucs.append(roc_auc_rank(y_perm, scores))
        aucs = np.asarray(aucs, dtype=float)
        auc_null = aucs
        p_perm = (1.0 + float((aucs >= auc).sum())) / (1.0 + len(aucs))

    # Bootstrap CIs (optional)
    boot = None
    auc_ci = (float("nan"), float("nan"))
    ap_ci = (float("nan"), float("nan"))
    if args.bootstrap and args.bootstrap > 0:
        boot = bootstrap_summary(
            y=y,
            scores=scores,
            rule=args.rule,
            min_sens=args.min_sensitivity,
            max_fpr=args.max_fpr,
            n_boot=args.bootstrap,
            seed=args.seed,
        )
        auc_ci = boot["auc_ci95"]
        ap_ci = boot["ap_ci95"]

    # ----------------------------
    # Write outputs (1D)
    # ----------------------------
    df.to_csv(outdir / "per_sample_alignment_metrics_with_truth.tsv", sep="\t", index=False)
    perf.to_csv(outdir / "threshold_scan_1d.tsv", sep="\t", index=False)

    # Bootstrap outputs
    if boot is not None:
        boot_df = pd.DataFrame({
            "auc": boot["aucs"],
            "ap": boot["aps"],
            "selected_threshold": boot["ths"],
            "sensitivity": boot["sens"],
            "specificity": boot["spec"],
            "precision": boot["prec"],
            "f1": boot["f1"],
            "fpr": boot["fpr"],
        })
        boot_df.to_csv(outdir / "bootstrap_samples_1d.tsv", sep="\t", index=False)

    # ----------------------------
    # Plots (1D)
    # ----------------------------
    save_roc_plot(
        outdir / "roc_alignment.png",
        roc_pts=roc_pts,
        auc=auc,
        auc_ci=auc_ci,
        metric=args.metric,
        n_pos=n_pos,
        n_neg=n_neg,
        p_perm=p_perm,
    )

    save_pr_plot(
        outdir / "pr_alignment.png",
        recall=recall,
        precision=precision,
        ap=ap_score,
        ap_ci=ap_ci,
        metric=args.metric,
        prevalence=prevalence,
    )

    rule_desc = args.rule
    if args.rule == "minscore_sens_fpr":
        rule_desc += f" (min_sens={args.min_sensitivity}, max_fpr={args.max_fpr})"
    best_t = float(best["threshold_score"]) if best is not None else None

    save_metrics_vs_threshold_plot(
        outdir / "metrics_vs_threshold_1d.png",
        perf=perf,
        metric=args.metric,
        best_threshold=best_t,
        rule_desc=rule_desc,
    )

    if boot is not None:
        save_bootstrap_hist(
            outdir / "bootstrap_auc_hist.png",
            boot["aucs"],
            xlabel="AUC",
            title=f"Bootstrap distribution of AUC (n={len(boot['aucs'])}, skipped={boot['skipped']})",
            vline=auc,
        )
        save_bootstrap_hist(
            outdir / "bootstrap_threshold_hist.png",
            boot["ths"],
            xlabel="Selected threshold (score)",
            title=f"Bootstrap distribution of selected threshold (rule={rule_desc})",
            vline=best_t if best_t is not None else float("nan"),
        )
        save_bootstrap_hist(
            outdir / "bootstrap_sensitivity_hist.png",
            boot["sens"],
            xlabel="Sensitivity at selected threshold",
            title="Bootstrap distribution of sensitivity at selected threshold",
            vline=float(best["sensitivity"]) if best is not None else float("nan"),
        )
        save_bootstrap_hist(
            outdir / "bootstrap_specificity_hist.png",
            boot["spec"],
            xlabel="Specificity at selected threshold",
            title="Bootstrap distribution of specificity at selected threshold",
            vline=float(best["specificity"]) if best is not None else float("nan"),
        )

    if auc_null is not None:
        save_bootstrap_hist(
            outdir / "auc_permutation_null.png",
            auc_null,
            xlabel="AUC under permuted labels",
            title=f"Permutation null of AUC (obs={auc:.3f}, p={p_perm:.3g}, n={len(auc_null)})",
            vline=auc,
        )

    # ----------------------------
    # 2D threshold search (optional)
    # ----------------------------
    best2d = None
    scan2d = None
    second_name = None

    if args.two_dim != "none":
        mapped = df_roc["cmv_mapped_reads_mapq"].to_numpy(dtype=float)

        if args.two_dim == "mapped_breadth":
            second = df_roc["cmv_breadth_mapq"].to_numpy(dtype=float)
            second_name = "breadth"
        elif args.two_dim == "mapped_covbases":
            second = df_roc["cmv_cov_bases_mapq"].to_numpy(dtype=float)
            second_name = "cov_bases"
        else:
            raise ValueError("Unexpected --two-dim value")

        scan2d, t_m, t_s = threshold_scan_2d(y, mapped, second, second_name=second_name)

        # rank + write top N (by chosen objective)
        if args.rule2d == "youden":
            obj = scan2d["youden_J"].fillna(-1e18)
            scan2d["objective"] = obj
        elif args.rule2d == "max_f1":
            obj = scan2d["F1"].fillna(-1e18)
            scan2d["objective"] = obj
        else:
            # for constraint rule, use objective=Youden for sorting, but selection is lexicographic within feasible set
            scan2d["objective"] = scan2d["youden_J"].fillna(-1e18)

        scan2d_sorted = scan2d.sort_values("objective", ascending=False).reset_index(drop=True)
        scan2d_sorted.head(int(args.top2d)).to_csv(outdir / f"threshold_scan_2d_top{args.top2d}_{args.two_dim}.tsv", sep="\t", index=False)
        scan2d.to_csv(outdir / f"threshold_scan_2d_full_{args.two_dim}.tsv", sep="\t", index=False)

        best2d = select_threshold_2d(
            scan2d=scan2d,
            rule=args.rule2d,
            min_sens=args.min_sensitivity,
            max_fpr=args.max_fpr,
            second_name=second_name,
        )

        # Heatmaps (Youden + F1 are usually most interpretable)
        extra = ""
        if args.rule2d == "minscore_sens_fpr":
            extra = f"\nConstraint: sens≥{args.min_sensitivity}, FPR≤{args.max_fpr}"
        save_2d_heatmap(
            outdir / f"heatmap_2d_youden_{args.two_dim}.png",
            scan2d=scan2d,
            t_m=t_m,
            t_s=t_s,
            second_name=second_name,
            value_col="youden_J",
            best_pair=best2d,
            title_extra=extra,
        )
        save_2d_heatmap(
            outdir / f"heatmap_2d_f1_{args.two_dim}.png",
            scan2d=scan2d,
            t_m=t_m,
            t_s=t_s,
            second_name=second_name,
            value_col="F1",
            best_pair=best2d,
            title_extra=extra,
        )

    # ----------------------------
    # report
    # ----------------------------
    with open(outdir / "report_alignment_thresholding.txt", "w") as f:
        f.write("Alignment-based CMV detection threshold optimisation\n\n")
        f.write(f"Truth samples with BAMs: n={len(y)} (pos={n_pos}, neg={n_neg}, prevalence={prevalence:.4g})\n\n")

        f.write("=== Per-BAM metric settings ===\n")
        f.write(f"MAPQ cutoff: {args.mapq}\n")
        f.write(f"Breadth method: {args.breadth_method} (depth>= {args.breadth_depth_min})\n\n")

        f.write("=== 1D Score evaluation ===\n")
        f.write(f"Score metric: {args.metric}\n")
        f.write(f"ROC AUC: {auc:.6g}\n")
        f.write(f"PR Average Precision (AP): {ap_score:.6g}\n")

        if p_perm is not None:
            f.write(f"Permutation p-value (one-sided, AUC>=obs): {p_perm:.6g} (n={args.permute})\n")

        if boot is not None:
            f.write("\n--- Bootstrap 95% CIs (percentile) ---\n")
            f.write(f"AUC 95% CI: {boot['auc_ci95'][0]:.6g} – {boot['auc_ci95'][1]:.6g}\n")
            f.write(f"AP  95% CI: {boot['ap_ci95'][0]:.6g} – {boot['ap_ci95'][1]:.6g}\n")
            f.write(f"Selected threshold 95% CI: {boot['th_ci95'][0]:.6g} – {boot['th_ci95'][1]:.6g}\n")
            f.write(f"Sensitivity 95% CI: {boot['sens_ci95'][0]:.6g} – {boot['sens_ci95'][1]:.6g}\n")
            f.write(f"Specificity 95% CI: {boot['spec_ci95'][0]:.6g} – {boot['spec_ci95'][1]:.6g}\n")
            f.write(f"Precision   95% CI: {boot['prec_ci95'][0]:.6g} – {boot['prec_ci95'][1]:.6g}\n")
            f.write(f"F1          95% CI: {boot['f1_ci95'][0]:.6g} – {boot['f1_ci95'][1]:.6g}\n")
            f.write(f"FPR         95% CI: {boot['fpr_ci95'][0]:.6g} – {boot['fpr_ci95'][1]:.6g}\n")
            f.write(f"Bootstrap iterations skipped (all-pos/all-neg resamples): {boot['skipped']}\n")

        f.write("\n=== 1D Threshold selection ===\n")
        if args.rule == "youden":
            f.write("Rule: maximise Youden's J (sens + spec - 1)\n")
        else:
            f.write(f"Rule: smallest threshold meeting constraints: sens≥{args.min_sensitivity}, FPR≤{args.max_fpr}\n")

        if best is None:
            f.write("Selected threshold: NONE (no feasible threshold found)\n")
        else:
            f.write(f"Selected threshold (score): {best['threshold_score']:.6g}\n")
            f.write(f"TP={best['TP']} TN={best['TN']} FP={best['FP']} FN={best['FN']}\n")
            f.write(f"Sensitivity={best['sensitivity']:.4g} Specificity={best['specificity']:.4g} FPR={best['FPR']:.4g}\n")
            f.write(f"Precision={best['precision']:.4g} F1={best['F1']:.4g} YoudenJ={best['youden_J']:.4g}\n")

        if args.two_dim != "none":
            f.write("\n=== 2D Threshold search ===\n")
            f.write(f"Mode: {args.two_dim} (positive if mapped_hq>=Tm AND {second_name}>=Ts)\n")
            f.write(f"Selection rule (2D): {args.rule2d}\n")
            if args.rule2d == "minscore_sens_fpr":
                f.write(f"Constraints: sens≥{args.min_sensitivity}, FPR≤{args.max_fpr}\n")
            if best2d is None:
                f.write("Selected 2D threshold pair: NONE\n")
            else:
                f.write(f"Selected Tm (mapped_hq): {best2d['Tm_mapped_hq']:.6g}\n")
                f.write(f"Selected Ts ({second_name}): {best2d[f'Ts_{second_name}']:.6g}\n")
                f.write(f"TP={best2d['TP']} TN={best2d['TN']} FP={best2d['FP']} FN={best2d['FN']}\n")
                f.write(f"Sensitivity={best2d['sensitivity']:.4g} Specificity={best2d['specificity']:.4g} FPR={best2d['FPR']:.4g}\n")
                f.write(f"Precision={best2d['precision']:.4g} F1={best2d['F1']:.4g} YoudenJ={best2d['youden_J']:.4g}\n")

    # Console summary
    print("✅ Done")
    print(f"Outputs in: {outdir}")
    print("Key outputs:")
    print("  - per_sample_alignment_metrics_with_truth.tsv")
    print("  - threshold_scan_1d.tsv")
    print("  - roc_alignment.png")
    print("  - pr_alignment.png")
    print("  - metrics_vs_threshold_1d.png")
    print("  - report_alignment_thresholding.txt")
    if boot is not None:
        print("  - bootstrap_samples_1d.tsv")
        print("  - bootstrap_auc_hist.png")
        print("  - bootstrap_threshold_hist.png")
        print("  - bootstrap_sensitivity_hist.png")
        print("  - bootstrap_specificity_hist.png")
    if auc_null is not None:
        print("  - auc_permutation_null.png")
    if args.two_dim != "none":
        print(f"  - threshold_scan_2d_full_{args.two_dim}.tsv")
        print(f"  - threshold_scan_2d_top{args.top2d}_{args.two_dim}.tsv")
        print(f"  - heatmap_2d_youden_{args.two_dim}.png")
        print(f"  - heatmap_2d_f1_{args.two_dim}.png")


if __name__ == "__main__":
    main()
