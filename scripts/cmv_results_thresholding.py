#!/usr/bin/env python3
"""
cmv_results_thresholding_v3_kofn.py

Adds comparison of 1/3, 2/3, and 3/3 multi-metric CMV detection rules.

Expected metric columns in --metrics-tsv:
  - log10_rpm_hq
  - cmv_breadth_mapq
  - blast_cmv_identity_top5_mean

Outputs:
  - merged_metrics_with_truth.tsv
  - threshold_table.tsv                  (for single score mode if used)
  - multi_metric_rule_comparison.tsv     (1/3 vs 2/3 vs 3/3)
  - report_results.txt

Figures:
  - fig1_score_distributions.png         (single-score mode)
  - fig2_roc_curve_step.png              (single-score mode)
  - fig3_pr_curve_step.png               (single-score mode)
  - fig4_metrics_vs_threshold.png        (single-score mode)
  - fig5_confusion_matrix.png            (single-score mode)
  - fig6_group_outcomes.png              (single-score mode)

  New:
  - fig12_confusion_matrix_1of3.png
  - fig13_confusion_matrix_2of3.png
  - fig14_confusion_matrix_3of3.png
  - fig15_kofn_rule_comparison.png
"""

import argparse
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------
# sample name canonicalisation
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


# ----------------------------
# core classification helpers
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
    recall = safe_div(tp, tp + fn)
    spec = safe_div(tn, tn + fp)
    prec = safe_div(tp, tp + fp)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    fpr = safe_div(fp, fp + tn)
    acc = safe_div(tp + tn, tp + tn + fp + fn)
    youden = recall + spec - 1 if (np.isfinite(recall) and np.isfinite(spec)) else float("nan")
    return recall, spec, prec, f1, fpr, acc, youden


# ----------------------------
# ROC / PR
# ----------------------------
def roc_curve_step(y_true, scores):
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores).astype(float)

    thr = np.unique(s)
    thr = np.sort(thr)
    thresholds = np.concatenate(([thr[-1] + 1e-12], thr[::-1], [thr[0] - 1e-12]))

    tpr = []
    fpr = []
    th_out = []

    for t in thresholds:
        yhat = (s >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y, yhat)
        rec, sp, pr, f1, fp_rate, acc, yj = metrics_from_counts(tp, tn, fp, fn)
        tpr.append(rec)
        fpr.append(fp_rate)
        th_out.append(float(t))

    fpr = np.asarray(fpr, float)
    tpr = np.asarray(tpr, float)
    th_out = np.asarray(th_out, float)
    order = np.argsort(fpr)
    return fpr[order], tpr[order], th_out[order]


def auc_trapz(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < 2:
        return float("nan")
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    return float(np.trapz(y, x))


def pr_curve_step(y_true, scores):
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores).astype(float)

    thr = np.unique(s)
    thr = np.sort(thr)
    thresholds = np.concatenate(([thr[-1] + 1e-12], thr[::-1], [thr[0] - 1e-12]))

    recs = []
    precs = []
    th_out = []

    for t in thresholds:
        yhat = (s >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y, yhat)
        rec, sp, pr, f1, fp_rate, acc, yj = metrics_from_counts(tp, tn, fp, fn)
        recs.append(rec)
        precs.append(pr)
        th_out.append(float(t))

    recs = np.asarray(recs, float)
    precs = np.asarray(precs, float)
    th_out = np.asarray(th_out, float)

    order = np.argsort(recs)
    return recs[order], precs[order], th_out[order]


def average_precision(y_true, scores):
    rec, prec, _ = pr_curve_step(y_true, scores)
    m = np.isfinite(rec) & np.isfinite(prec)
    rec = rec[m]
    prec = prec[m]
    if len(rec) < 2:
        return float("nan")
    order = np.argsort(rec)
    rec = rec[order]
    prec = prec[order]
    return float(np.trapz(prec, rec))


# ----------------------------
# Threshold table + selection rules
# ----------------------------
def threshold_table(y_true, scores):
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores).astype(float)

    thr = np.unique(s)
    thr = np.sort(thr)
    thresholds = np.concatenate(([thr[-1] + 1e-12], thr[::-1], [thr[0] - 1e-12]))

    rows = []
    for t in thresholds:
        yhat = (s >= t).astype(int)
        tp, tn, fp, fn = confusion_counts(y, yhat)
        rec, sp, pr, f1, fp_rate, acc, yj = metrics_from_counts(tp, tn, fp, fn)
        rows.append({
            "threshold": float(t),
            "TP": tp, "TN": tn, "FP": fp, "FN": fn,
            "recall": rec,
            "specificity": sp,
            "precision": pr,
            "F1": f1,
            "FPR": fp_rate,
            "accuracy": acc,
            "youden_J": yj,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("threshold", ascending=False).reset_index(drop=True)
    return df


def pick_threshold(df: pd.DataFrame, rule: str,
                   min_recall: float, min_precision: float,
                   max_fpr: float):
    d = df.copy()

    if rule == "youden":
        d["obj"] = d["youden_J"].fillna(-1e18)
        return d.loc[d["obj"].idxmax()].to_dict()

    if rule == "max_f1":
        d["obj"] = d["F1"].fillna(-1e18)
        return d.loc[d["obj"].idxmax()].to_dict()

    if rule == "max_precision_at_min_recall":
        feas = d[d["recall"] >= min_recall].copy()
        if feas.empty:
            return None
        feas["obj"] = feas["precision"].fillna(-1e18)
        feas = feas.sort_values(["obj", "recall", "threshold"], ascending=[False, False, True])
        return feas.iloc[0].to_dict()

    if rule == "max_recall_at_min_precision":
        feas = d[d["precision"] >= min_precision].copy()
        if feas.empty:
            return None
        feas["obj"] = feas["recall"].fillna(-1e18)
        feas = feas.sort_values(["obj", "precision", "threshold"], ascending=[False, False, True])
        return feas.iloc[0].to_dict()

    if rule == "min_threshold_meeting_constraints":
        feas = d[(d["recall"] >= min_recall) & (d["FPR"] <= max_fpr)].copy()
        if feas.empty:
            return None
        feas = feas.sort_values("threshold", ascending=True)
        return feas.iloc[0].to_dict()

    raise ValueError(f"Unknown rule: {rule}")


# ----------------------------
# Fisher exact test
# ----------------------------
def _hypergeom_p(a, b, c, d):
    row1 = a + b
    row2 = c + d
    col1 = a + c
    col2 = b + d
    n = row1 + row2
    num = math.comb(col1, a) * math.comb(col2, row1 - a)
    den = math.comb(n, row1)
    return num / den


def fishers_exact_2x2(a, b, c, d, alternative="two-sided"):
    or_num = (a + 0.5) * (d + 0.5)
    or_den = (b + 0.5) * (c + 0.5)
    odds_ratio = or_num / or_den

    row1 = a + b
    col1 = a + c
    col2 = b + d
    n = a + b + c + d

    a_min = max(0, row1 - col2)
    a_max = min(row1, col1)
    p_obs = _hypergeom_p(a, b, c, d)

    ps = []
    for a2 in range(a_min, a_max + 1):
        b2 = row1 - a2
        c2 = col1 - a2
        d2 = n - a2 - b2 - c2
        p = _hypergeom_p(a2, b2, c2, d2)
        ps.append((a2, p))

    if alternative == "two-sided":
        p_val = sum(p for _, p in ps if p <= p_obs + 1e-15)
    elif alternative == "greater":
        p_val = sum(p for a2, p in ps if a2 >= a)
    elif alternative == "less":
        p_val = sum(p for a2, p in ps if a2 <= a)
    else:
        raise ValueError("alternative must be two-sided/greater/less")

    p_val = min(1.0, max(0.0, p_val))
    return float(odds_ratio), float(p_val)


# ----------------------------
# bootstrap + permutation
# ----------------------------
def percentile_ci(x, alpha=0.05):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (float("nan"), float("nan"))
    lo = np.percentile(x, 100 * (alpha / 2))
    hi = np.percentile(x, 100 * (1 - alpha / 2))
    return float(lo), float(hi)


def bootstrap(y, scores, rule, min_recall, min_precision, max_fpr, n_boot=2000, seed=1):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, int)
    scores = np.asarray(scores, float)
    n = len(y)

    aucs, aps, ths, recs, precs, f1s, specs = [], [], [], [], [], [], []
    skipped = 0

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        sb = scores[idx]
        if (yb == 1).sum() == 0 or (yb == 0).sum() == 0:
            skipped += 1
            continue

        fpr, tpr, _ = roc_curve_step(yb, sb)
        aucs.append(auc_trapz(fpr, tpr))
        aps.append(average_precision(yb, sb))

        tab = threshold_table(yb, sb)
        best = pick_threshold(tab, rule, min_recall, min_precision, max_fpr)
        if best is None:
            ths.append(float("nan"))
            recs.append(float("nan"))
            precs.append(float("nan"))
            f1s.append(float("nan"))
            specs.append(float("nan"))
        else:
            ths.append(best["threshold"])
            recs.append(best["recall"])
            precs.append(best["precision"])
            f1s.append(best["F1"])
            specs.append(best["specificity"])

    out = {
        "aucs": np.asarray(aucs, float),
        "aps": np.asarray(aps, float),
        "ths": np.asarray(ths, float),
        "recall": np.asarray(recs, float),
        "precision": np.asarray(precs, float),
        "F1": np.asarray(f1s, float),
        "specificity": np.asarray(specs, float),
        "skipped": int(skipped),
    }
    out["auc_ci95"] = percentile_ci(out["aucs"])
    out["ap_ci95"] = percentile_ci(out["aps"])
    out["th_ci95"] = percentile_ci(out["ths"])
    out["recall_ci95"] = percentile_ci(out["recall"])
    out["precision_ci95"] = percentile_ci(out["precision"])
    out["f1_ci95"] = percentile_ci(out["F1"])
    out["spec_ci95"] = percentile_ci(out["specificity"])
    return out


def permutation_test(y, scores, n_perm=2000, seed=1):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, int)
    scores = np.asarray(scores, float)

    fpr, tpr, _ = roc_curve_step(y, scores)
    auc_obs = auc_trapz(fpr, tpr)
    ap_obs = average_precision(y, scores)

    auc_null, ap_null = [], []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        fprp, tprp, _ = roc_curve_step(yp, scores)
        auc_null.append(auc_trapz(fprp, tprp))
        ap_null.append(average_precision(yp, scores))

    auc_null = np.asarray(auc_null, float)
    ap_null = np.asarray(ap_null, float)

    p_auc = (1.0 + float((auc_null >= auc_obs).sum())) / (1.0 + len(auc_null))
    p_ap = (1.0 + float((ap_null >= ap_obs).sum())) / (1.0 + len(ap_null))

    return auc_obs, ap_obs, auc_null, ap_null, float(p_auc), float(p_ap)


# ----------------------------
# samtools helpers (optional BAM metrics)
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


def breadth_of_coverage(bam: Path, mapq: int, depth_min: int = 1):
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
# IO
# ----------------------------
def load_truth(truth_path: Path) -> pd.DataFrame:
    df = pd.read_csv(truth_path, sep="\t")
    if "sample" not in df.columns:
        raise ValueError("Truth TSV must contain 'sample'.")
    df["base"] = df["sample"].astype(str).apply(canonical_sample)

    if "truth_cmv_present" in df.columns:
        y = df["truth_cmv_present"]
    elif "cmv_present" in df.columns:
        y = df["cmv_present"]
    else:
        raise ValueError("Truth TSV must contain 'truth_cmv_present' or 'cmv_present'.")

    df["truth_cmv_present"] = pd.to_numeric(y, errors="coerce").fillna(0).astype(int)

    if "group" not in df.columns:
        df["group"] = "UNKNOWN"

    df = df.sort_values("base").drop_duplicates("base", keep="first").reset_index(drop=True)
    return df


def find_bams(cmv_bam_dir: Path):
    bams = sorted([p for p in cmv_bam_dir.glob("*.bam") if not str(p).endswith(".bai")])
    if not bams:
        raise FileNotFoundError(f"No .bam files found in: {cmv_bam_dir}")
    return bams


# ----------------------------
# plotting
# ----------------------------
def save_hist(path: Path, values: np.ndarray, xlabel: str, title: str, vline=None):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    plt.figure(figsize=(6.4, 5.2))
    if v.size:
        plt.hist(v, bins=45)
    if vline is not None and np.isfinite(vline):
        plt.axvline(vline, linestyle="--", label="Observed")
        plt.legend(frameon=True)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_distributions(path: Path, df: pd.DataFrame, score_col: str):
    score = df[score_col].astype(float).to_numpy()
    score_t = np.log10(1 + np.clip(score, 0, None))
    y = df["truth_cmv_present"].astype(int).to_numpy()
    groups = df["group"].astype(str).to_numpy()

    preferred = [
        "NEG_HUMAN_ONLY", "NEG_OTHER_VIRUS", "NEG_CONTAM",
        "POS_ROC_LADDER", "POS_CMV_PLUS_OTHER", "POS_MIXED_STRAINS", "POS_MUT_SPIKE"
    ]
    uniq = list(dict.fromkeys(groups))
    ordered = [g for g in preferred if g in uniq] + [g for g in uniq if g not in preferred]

    rng = np.random.default_rng(1)

    plt.figure(figsize=(11.5, 5.3))

    ax1 = plt.subplot(1, 2, 1)
    jitter = rng.normal(0, 0.04, size=len(score_t))
    ax1.scatter(y + jitter, score_t, s=18)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["Truth negative", "Truth positive"])
    ax1.set_xlabel("Truth label")
    ax1.set_ylabel(f"log10(1 + {score_col})")
    ax1.set_title("Score distribution by truth label")
    ax1.text(0.02, 0.98, f"n_neg={(y==0).sum()}\nn_pos={(y==1).sum()}",
             transform=ax1.transAxes, va="top")

    ax2 = plt.subplot(1, 2, 2)
    g_to_x = {g: i for i, g in enumerate(ordered)}
    xs = np.array([g_to_x[g] for g in groups], float)
    jitter2 = rng.normal(0, 0.08, size=len(xs))
    ax2.scatter(xs + jitter2, score_t, s=18)
    ax2.set_xticks(range(len(ordered)))
    ax2.set_xticklabels(ordered, rotation=45, ha="right")
    ax2.set_xlabel("Synthetic sample group")
    ax2.set_ylabel(f"log10(1 + {score_col})")
    ax2.set_title("Score distribution by synthetic group")

    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_roc(path: Path, fpr, tpr, auc, auc_ci=None, p=None, metric_name="score"):
    plt.figure(figsize=(6.4, 5.2))
    plt.step(fpr, tpr, where="post", label="ROC (step)")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Chance")
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    title = f"ROC (metric={metric_name}) | AUC={auc:.3f}"
    if auc_ci is not None and np.isfinite(auc_ci[0]) and np.isfinite(auc_ci[1]):
        title += f" (95% CI {auc_ci[0]:.3f}–{auc_ci[1]:.3f})"
    if p is not None and np.isfinite(p):
        title += f" | permutation p={p:.3g}"
    plt.title(title)
    plt.legend(loc="lower right", frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_pr(path: Path, recall, precision, ap, ap_ci=None, p=None, prevalence=None, metric_name="score"):
    plt.figure(figsize=(6.4, 5.2))
    plt.step(recall, precision, where="post", label="PR (step)")
    if prevalence is not None and np.isfinite(prevalence):
        plt.axhline(prevalence, linestyle="--", label=f"Baseline prevalence={prevalence:.3f}")
    plt.xlabel("Recall (Sensitivity)")
    plt.ylabel("Precision (PPV)")
    plt.ylim(0, 1.02)
    plt.xlim(0, 1.0)
    title = f"PR (metric={metric_name}) | AP={ap:.3f}"
    if ap_ci is not None and np.isfinite(ap_ci[0]) and np.isfinite(ap_ci[1]):
        title += f" (95% CI {ap_ci[0]:.3f}–{ap_ci[1]:.3f})"
    if p is not None and np.isfinite(p):
        title += f" | permutation p={p:.3g}"
    plt.title(title)
    plt.legend(loc="lower left", frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_metrics_vs_threshold(path: Path, tab: pd.DataFrame, chosen_t=None, rule_desc=""):
    plt.figure(figsize=(8.0, 5.2))
    d = tab.sort_values("threshold", ascending=True)
    plt.plot(d["threshold"], d["recall"], label="Recall")
    plt.plot(d["threshold"], d["specificity"], label="Specificity")
    plt.plot(d["threshold"], d["precision"], label="Precision")
    plt.plot(d["threshold"], d["F1"], label="F1")
    plt.plot(d["threshold"], d["youden_J"], label="Youden J")

    if chosen_t is not None and np.isfinite(chosen_t):
        plt.axvline(chosen_t, linestyle="--", label=f"Selected T={chosen_t:.4g}")

    plt.xlabel("Threshold T")
    plt.ylabel("Metric value")
    title = "Metrics vs threshold"
    if rule_desc:
        title += f"\nSelection rule: {rule_desc}"
    plt.title(title)
    plt.legend(frameon=True, ncol=2)
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_confusion(path: Path, tp, tn, fp, fn, title):
    mat = np.array([[tn, fp],
                    [fn, tp]], float)
    plt.figure(figsize=(5.6, 4.8))
    im = plt.imshow(mat, aspect="auto")
    plt.colorbar(im, label="Number of samples")
    plt.xticks([0, 1], ["Pred NEG", "Pred POS"])
    plt.yticks([0, 1], ["True NEG", "True POS"])
    plt.title(title)
    for (i, j), v in np.ndenumerate(mat):
        plt.text(j, i, f"{int(v)}", ha="center", va="center")
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_group_outcomes(path: Path, df: pd.DataFrame, score_col: str, threshold: float):
    d = df.copy()
    y = d["truth_cmv_present"].astype(int).to_numpy()
    s = d[score_col].astype(float).to_numpy()
    yhat = (s >= threshold).astype(int)

    outcome = []
    for yt, yp in zip(y, yhat):
        if yt == 1 and yp == 1:
            outcome.append("TP")
        elif yt == 0 and yp == 0:
            outcome.append("TN")
        elif yt == 0 and yp == 1:
            outcome.append("FP")
        else:
            outcome.append("FN")
    d["outcome"] = outcome

    counts = d.groupby(["group", "outcome"]).size().unstack(fill_value=0)

    for col in ["TP", "FP", "FN", "TN"]:
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[["TP", "FP", "FN", "TN"]]

    plt.figure(figsize=(10.5, 5.2))
    x = np.arange(len(counts.index))
    bottom = np.zeros(len(counts.index))
    for col in ["TP", "FP", "FN", "TN"]:
        vals = counts[col].to_numpy()
        plt.bar(x, vals, bottom=bottom, label=col)
        bottom += vals

    plt.xticks(x, counts.index, rotation=45, ha="right")
    plt.ylabel("Number of samples")
    plt.title(f"Per-group outcomes at threshold T={threshold:.4g}")
    plt.legend(frameon=True, ncol=4)
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


def plot_kofn_comparison(path: Path, comp_df: pd.DataFrame):
    d = comp_df.copy()
    rules = d["rule"].tolist()

    plt.figure(figsize=(8.5, 5.2))
    x = np.arange(len(rules))
    width = 0.18

    plt.bar(x - 1.5 * width, d["recall"], width=width, label="Recall")
    plt.bar(x - 0.5 * width, d["specificity"], width=width, label="Specificity")
    plt.bar(x + 0.5 * width, d["precision"], width=width, label="Precision")
    plt.bar(x + 1.5 * width, d["F1"], width=width, label="F1")

    plt.xticks(x, rules)
    plt.ylim(0, 1.05)
    plt.ylabel("Metric value")
    plt.title("Comparison of 1/3, 2/3, and 3/3 CMV detection rules")
    plt.legend(frameon=True, ncol=2)
    plt.tight_layout()
    plt.savefig(path, dpi=240)
    plt.close()


# ----------------------------
# multi-metric k-of-n comparison
# ----------------------------
def evaluate_kofn_rule(df: pd.DataFrame,
                       rpm_col: str,
                       breadth_col: str,
                       blast_col: str,
                       rpm_thr: float,
                       breadth_thr: float,
                       blast_thr: float,
                       k: int) -> dict:
    d = df.copy()

    d["pass_rpm"] = pd.to_numeric(d[rpm_col], errors="coerce") >= rpm_thr
    d["pass_breadth"] = pd.to_numeric(d[breadth_col], errors="coerce") >= breadth_thr
    d["pass_blast"] = pd.to_numeric(d[blast_col], errors="coerce") >= blast_thr

    d["n_pass"] = d[["pass_rpm", "pass_breadth", "pass_blast"]].sum(axis=1)
    d["pred_kofn"] = (d["n_pass"] >= k).astype(int)

    y_true = d["truth_cmv_present"].astype(int).to_numpy()
    y_pred = d["pred_kofn"].astype(int).to_numpy()

    tp, tn, fp, fn = confusion_counts(y_true, y_pred)
    recall, spec, prec, f1, fpr, acc, youden = metrics_from_counts(tp, tn, fp, fn)
    or_val, p_fisher = fishers_exact_2x2(a=tp, b=fn, c=fp, d=tn, alternative="two-sided")

    return {
        "rule": f"{k}/3",
        "k": k,
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "recall": recall,
        "specificity": spec,
        "precision": prec,
        "F1": f1,
        "FPR": fpr,
        "accuracy": acc,
        "youden_J": youden,
        "odds_ratio": or_val,
        "fisher_p": p_fisher,
        "df_with_preds": d,
    }


# ----------------------------
# main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="metadata.tsv from synthetic generator")
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--metrics-tsv", default="", help="Precomputed metrics TSV")
    ap.add_argument("--cmv-bam-dir", default="", help="Dir of CMV BAMs (optional)")
    ap.add_argument("--mapq", type=int, default=30)
    ap.add_argument("--breadth-depth-min", type=int, default=1)

    ap.add_argument("--score", choices=["mapped_hq", "breadth", "cov_bases", "combo"], default="mapped_hq")
    ap.add_argument("--rule",
                    choices=[
                        "youden",
                        "max_f1",
                        "max_precision_at_min_recall",
                        "max_recall_at_min_precision",
                        "min_threshold_meeting_constraints",
                    ],
                    default="max_f1")
    ap.add_argument("--min-recall", type=float, default=0.99)
    ap.add_argument("--min-precision", type=float, default=0.99)
    ap.add_argument("--max-fpr", type=float, default=0.05)

    ap.add_argument("--permute", type=int, default=2000)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)

    # new multi-metric arguments
    ap.add_argument("--rpm-col", default="log10_rpm_hq")
    ap.add_argument("--breadth-col", default="cmv_breadth_mapq")
    ap.add_argument("--blast-col", default="blast_cmv_identity_top5_mean")

    ap.add_argument("--rpm-threshold", type=float, default=2.0)
    ap.add_argument("--breadth-threshold", type=float, default=0.0061)
    ap.add_argument("--blast-threshold", type=float, default=87.0)

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    truth = load_truth(Path(args.truth))

    # ----------------------------
    # Load metrics
    # ----------------------------
    if args.metrics_tsv:
        mdf = pd.read_csv(args.metrics_tsv, sep="\t")
        if "base" not in mdf.columns:
            if "sample" not in mdf.columns:
                raise ValueError("metrics TSV must have base or sample column.")
            mdf["base"] = mdf["sample"].astype(str).apply(canonical_sample)
    else:
        if not args.cmv_bam_dir:
            raise ValueError("Provide --metrics-tsv or --cmv-bam-dir.")
        bams = find_bams(Path(args.cmv_bam_dir))
        rows = []
        for bam in bams:
            sample = re.sub(r"\.bam$", "", bam.name, flags=re.IGNORECASE)
            base = canonical_sample(sample)

            mapped_hq = samtools_count(bam, mapped_only=True, mapq=args.mapq)
            breadth, cov_bases, ref_bases = breadth_of_coverage(
                bam, mapq=args.mapq, depth_min=args.breadth_depth_min
            )

            rows.append({
                "base": base,
                "sample": sample,
                "cmv_mapped_reads_mapq": int(mapped_hq),
                "cmv_breadth_mapq": float(breadth),
                "cmv_cov_bases_mapq": int(cov_bases),
                "cmv_ref_bases": int(ref_bases),
                "bam": str(bam),
            })
        mdf = pd.DataFrame(rows)

    df = truth.merge(mdf, on="base", how="left")
    df.to_csv(outdir / "merged_metrics_with_truth.tsv", sep="\t", index=False)

    # ----------------------------
    # Single-score section (kept from your original script)
    # ----------------------------
    if args.score == "mapped_hq":
        score_col = "cmv_mapped_reads_mapq"
        df["score_used"] = pd.to_numeric(df[score_col], errors="coerce")
    elif args.score == "breadth":
        score_col = "cmv_breadth_mapq"
        df["score_used"] = pd.to_numeric(df[score_col], errors="coerce")
    elif args.score == "cov_bases":
        score_col = "cmv_cov_bases_mapq"
        df["score_used"] = pd.to_numeric(df[score_col], errors="coerce")
    else:
        df["score_used"] = np.log1p(pd.to_numeric(df["cmv_mapped_reads_mapq"], errors="coerce").fillna(0)) + \
                           10.0 * pd.to_numeric(df["cmv_breadth_mapq"], errors="coerce").fillna(0)
        score_col = "score_used"

    df_ana = df.dropna(subset=["truth_cmv_present", "score_used"]).copy()
    y = df_ana["truth_cmv_present"].astype(int).to_numpy()
    scores = df_ana["score_used"].astype(float).to_numpy()

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    prevalence = safe_div(n_pos, len(y))

    plot_distributions(outdir / "fig1_score_distributions.png", df_ana, "score_used")

    fpr, tpr, _ = roc_curve_step(y, scores)
    auc_obs = auc_trapz(fpr, tpr)

    rec, prec, _ = pr_curve_step(y, scores)
    ap_obs = average_precision(y, scores)

    tab = threshold_table(y, scores)
    tab.to_csv(outdir / "threshold_table.tsv", sep="\t", index=False)

    best = pick_threshold(tab, args.rule, args.min_recall, args.min_precision, args.max_fpr)
    if best is None:
        raise RuntimeError("No feasible single-score threshold found.")

    T = float(best["threshold"])
    tp, tn, fp_, fn = int(best["TP"]), int(best["TN"]), int(best["FP"]), int(best["FN"])
    recall, spec, precision, f1, fpr_best, acc, youden = metrics_from_counts(tp, tn, fp_, fn)

    rule_desc = args.rule
    if args.rule in {"max_precision_at_min_recall", "min_threshold_meeting_constraints"}:
        rule_desc += f" (min_recall={args.min_recall})"
    if args.rule in {"max_recall_at_min_precision"}:
        rule_desc += f" (min_precision={args.min_precision})"
    if args.rule == "min_threshold_meeting_constraints":
        rule_desc += f" (max_fpr={args.max_fpr})"

    plot_metrics_vs_threshold(outdir / "fig4_metrics_vs_threshold.png", tab, chosen_t=T, rule_desc=rule_desc)
    plot_confusion(outdir / "fig5_confusion_matrix.png", tp, tn, fp_, fn,
                   title=f"Confusion matrix at selected threshold T={T:.4g}")
    plot_group_outcomes(outdir / "fig6_group_outcomes.png", df_ana, "score_used", threshold=T)

    odds_ratio, p_fisher = fishers_exact_2x2(a=tp, b=fn, c=fp_, d=tn, alternative="two-sided")

    p_auc = p_ap = None
    auc_null = ap_null = None
    if args.permute and args.permute > 0:
        auc_obs2, ap_obs2, auc_null, ap_null, p_auc, p_ap = permutation_test(
            y, scores, n_perm=args.permute, seed=args.seed
        )
        auc_obs = auc_obs2
        ap_obs = ap_obs2

        save_hist(outdir / "fig7_perm_null_auc.png", auc_null,
                  xlabel="AUC under permuted labels",
                  title=f"Permutation null (AUC) | observed={auc_obs:.3f}, p={p_auc:.3g}",
                  vline=auc_obs)

        save_hist(outdir / "fig8_perm_null_ap.png", ap_null,
                  xlabel="AP under permuted labels",
                  title=f"Permutation null (AP) | observed={ap_obs:.3f}, p={p_ap:.3g}",
                  vline=ap_obs)

    boot = None
    auc_ci = ap_ci = None
    if args.bootstrap and args.bootstrap > 0:
        boot = bootstrap(
            y, scores,
            rule=args.rule,
            min_recall=args.min_recall,
            min_precision=args.min_precision,
            max_fpr=args.max_fpr,
            n_boot=args.bootstrap,
            seed=args.seed
        )
        auc_ci = boot["auc_ci95"]
        ap_ci = boot["ap_ci95"]

        bdf = pd.DataFrame({
            "auc": boot["aucs"],
            "ap": boot["aps"],
            "threshold": boot["ths"],
            "recall": boot["recall"],
            "precision": boot["precision"],
            "F1": boot["F1"],
            "specificity": boot["specificity"],
        })
        bdf.to_csv(outdir / "bootstrap_samples.tsv", sep="\t", index=False)

        save_hist(outdir / "fig9_bootstrap_auc.png", boot["aucs"],
                  xlabel="AUC", title=f"Bootstrap AUC (n={len(boot['aucs'])})",
                  vline=auc_obs)

        save_hist(outdir / "fig10_bootstrap_ap.png", boot["aps"],
                  xlabel="AP", title=f"Bootstrap AP (n={len(boot['aps'])})",
                  vline=ap_obs)

        save_hist(outdir / "fig11_bootstrap_threshold.png", boot["ths"],
                  xlabel="Selected threshold", title=f"Bootstrap threshold | rule={args.rule}",
                  vline=T)

    plot_roc(outdir / "fig2_roc_curve_step.png", fpr, tpr, auc_obs, auc_ci=auc_ci, p=p_auc, metric_name=args.score)
    plot_pr(outdir / "fig3_pr_curve_step.png", rec, prec, ap_obs, ap_ci=ap_ci, p=p_ap,
            prevalence=prevalence, metric_name=args.score)

    # ----------------------------
    # New: multi-metric 1/3 vs 2/3 vs 3/3 comparison
    # ----------------------------
    needed_cols = [args.rpm_col, args.breadth_col, args.blast_col, "truth_cmv_present"]
    missing = [c for c in needed_cols if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns for k-of-n comparison: "
            + ", ".join(missing)
        )

    df_k = df.dropna(subset=needed_cols).copy()

    res_1 = evaluate_kofn_rule(
        df_k, args.rpm_col, args.breadth_col, args.blast_col,
        args.rpm_threshold, args.breadth_threshold, args.blast_threshold, k=1
    )
    res_2 = evaluate_kofn_rule(
        df_k, args.rpm_col, args.breadth_col, args.blast_col,
        args.rpm_threshold, args.breadth_threshold, args.blast_threshold, k=2
    )
    res_3 = evaluate_kofn_rule(
        df_k, args.rpm_col, args.breadth_col, args.blast_col,
        args.rpm_threshold, args.breadth_threshold, args.blast_threshold, k=3
    )

    comp_rows = []
    for res in [res_1, res_2, res_3]:
        comp_rows.append({
            "rule": res["rule"],
            "TP": res["TP"],
            "TN": res["TN"],
            "FP": res["FP"],
            "FN": res["FN"],
            "recall": res["recall"],
            "specificity": res["specificity"],
            "precision": res["precision"],
            "F1": res["F1"],
            "FPR": res["FPR"],
            "accuracy": res["accuracy"],
            "youden_J": res["youden_J"],
            "odds_ratio": res["odds_ratio"],
            "fisher_p": res["fisher_p"],
        })

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(outdir / "multi_metric_rule_comparison.tsv", sep="\t", index=False)

    plot_confusion(outdir / "fig12_confusion_matrix_1of3.png",
                   res_1["TP"], res_1["TN"], res_1["FP"], res_1["FN"],
                   title="Confusion matrix: 1 of 3 rule")

    plot_confusion(outdir / "fig13_confusion_matrix_2of3.png",
                   res_2["TP"], res_2["TN"], res_2["FP"], res_2["FN"],
                   title="Confusion matrix: 2 of 3 rule")

    plot_confusion(outdir / "fig14_confusion_matrix_3of3.png",
                   res_3["TP"], res_3["TN"], res_3["FP"], res_3["FN"],
                   title="Confusion matrix: 3 of 3 rule")

    plot_kofn_comparison(outdir / "fig15_kofn_rule_comparison.png", comp_df)

    # ----------------------------
    # Report
    # ----------------------------
    report_path = outdir / "report_results.txt"
    with open(report_path, "w") as f:
        f.write("Synthetic CMV dataset — Threshold optimisation results (v3)\n\n")
        f.write(f"n={len(y)} | pos={n_pos} | neg={n_neg} | prevalence={prevalence:.4g}\n\n")

        f.write("Single-score discrimination:\n")
        f.write(f"  ROC AUC = {auc_obs:.4g}\n")
        f.write(f"  PR  AP  = {ap_obs:.4g}\n")
        if auc_ci is not None:
            f.write(f"  AUC 95% CI: {auc_ci[0]:.4g}–{auc_ci[1]:.4g}\n")
        if ap_ci is not None:
            f.write(f"  AP  95% CI: {ap_ci[0]:.4g}–{ap_ci[1]:.4g}\n")
        if p_auc is not None:
            f.write(f"  Permutation p (AUC): {p_auc:.4g}\n")
        if p_ap is not None:
            f.write(f"  Permutation p (AP):  {p_ap:.4g}\n")
        f.write("\n")

        f.write("Single-score threshold selection:\n")
        f.write(f"  Rule: {rule_desc}\n")
        f.write(f"  Selected threshold T = {T:.6g}\n")
        f.write(f"  Confusion: TP={tp} TN={tn} FP={fp_} FN={fn}\n")
        f.write(f"  Recall      = {recall:.4g}\n")
        f.write(f"  Precision   = {precision:.4g}\n")
        f.write(f"  Specificity = {spec:.4g}\n")
        f.write(f"  F1          = {f1:.4g}\n")
        f.write(f"  Fisher p    = {p_fisher:.4g}\n\n")

        f.write("Multi-metric k-of-n comparison:\n")
        f.write(f"  Thresholds used: {args.rpm_col}>={args.rpm_threshold}, "
                f"{args.breadth_col}>={args.breadth_threshold}, "
                f"{args.blast_col}>={args.blast_threshold}\n\n")

        for _, row in comp_df.iterrows():
            f.write(f"  Rule {row['rule']}:\n")
            f.write(f"    TP={int(row['TP'])} TN={int(row['TN'])} FP={int(row['FP'])} FN={int(row['FN'])}\n")
            f.write(f"    Recall={row['recall']:.4g} Specificity={row['specificity']:.4g} "
                    f"Precision={row['precision']:.4g} F1={row['F1']:.4g}\n")
            f.write(f"    Fisher p={row['fisher_p']:.4g}\n\n")

        f.write("Key new outputs:\n")
        f.write("  fig12_confusion_matrix_1of3.png\n")
        f.write("  fig13_confusion_matrix_2of3.png\n")
        f.write("  fig14_confusion_matrix_3of3.png\n")
        f.write("  fig15_kofn_rule_comparison.png\n")
        f.write("  multi_metric_rule_comparison.tsv\n")

    print("✅ Done")
    print(f"Outputs in: {outdir}")
    print("k-of-n comparison:")
    print(comp_df.to_string(index=False))


if __name__ == "__main__":
    main()