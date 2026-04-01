#!/usr/bin/env python3
"""
c_validate_final_threshold_resampling.py

Bootstrap + permutation validation for a FINAL fixed CMV detection rule:
  CMV positive if >= k of 3 metrics pass:
    - log10 RPM (HQ CMV reads) >= thr_log10_rpm
    - CMV breadth (MAPQ-filtered) >= thr_breadth
    - BLAST mean % identity (top 5 hits) >= thr_blast

What this gives you (dissertation-friendly):
  - Bootstrap distributions of performance metrics + 95% CI
  - Permutation test p-value showing performance >> chance
  - Figures suitable for Results section

Inputs:
  --metrics metrics_table.tsv (from your thresholding pipeline)
  --thresholds-json thresholds.json (optional; if provided, uses thresholds from it)
  --outdir output folder

Outputs:
  Fig5_bootstrap_metric_distributions.png
  Fig6_permutation_null_F1.png
  resampling_summary.tsv

No sklearn/scipy required.
"""

import argparse
import json
from pathlib import Path

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
        f"Available columns (first 80): {list(df.columns)[:80]}"
    )


def safe_div(a, b):
    return a / b if b else 0.0


def confusion_counts(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn


def metrics_from_counts(tp, tn, fp, fn):
    recall = safe_div(tp, tp + fn)          # sensitivity
    spec = safe_div(tn, tn + fp)            # specificity
    prec = safe_div(tp, tp + fp)            # precision (PPV)
    f1 = safe_div(2 * prec * recall, prec + recall)
    bal_acc = 0.5 * (recall + spec)
    return recall, spec, prec, f1, bal_acc


def load_thresholds_json(path_json, defaults):
    thr_rpm, thr_breadth, thr_blast, k_required = defaults
    obj = json.loads(Path(path_json).read_text())

    # Support a few possible layouts
    if isinstance(obj, dict):
        # nested rule blocks
        for key in ["final_rule", "best_rule", "selected_rule", "rule"]:
            if key in obj and isinstance(obj[key], dict):
                sec = obj[key]
                thr_rpm = float(sec.get("thr_log10_rpm_hq", sec.get("thr_rpm", thr_rpm)))
                thr_breadth = float(sec.get("thr_cmv_breadth_mapq", sec.get("thr_breadth", thr_breadth)))
                thr_blast = float(sec.get("thr_blast_identity", sec.get("thr_blast", thr_blast)))
                k_required = int(sec.get("k_required", sec.get("k", k_required)))

        # flat keys
        thr_rpm = float(obj.get("thr_log10_rpm_hq", obj.get("thr_rpm", thr_rpm)))
        thr_breadth = float(obj.get("thr_cmv_breadth_mapq", obj.get("thr_breadth", thr_breadth)))
        thr_blast = float(obj.get("thr_blast_identity", obj.get("thr_blast", thr_blast)))
        k_required = int(obj.get("k_required", obj.get("k", k_required)))

    return thr_rpm, thr_breadth, thr_blast, k_required


def apply_rule(df, col_rpm, col_breadth, col_blast, thr_rpm, thr_breadth, thr_blast, k_required):
    """Apply k-of-3 rule; missing values count as FAIL (conservative, clinical-friendly)."""
    pass_sum = (
        (df[col_rpm] >= thr_rpm).fillna(False).astype(int)
        + (df[col_breadth] >= thr_breadth).fillna(False).astype(int)
        + (df[col_blast] >= thr_blast).fillna(False).astype(int)
    )
    return (pass_sum >= k_required).astype(int)


def percentile_ci(x, lo=2.5, hi=97.5):
    x = np.asarray(x, dtype=float)
    return float(np.percentile(x, lo)), float(np.percentile(x, hi))


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Bootstrap + permutation validation for the final fixed multi-metric CMV threshold."
    )
    ap.add_argument("--metrics", required=True, help="metrics_table.tsv")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--thresholds-json", default="", help="Optional thresholds.json")

    ap.add_argument("--thr-log10-rpm", type=float, default=2.0)
    ap.add_argument("--thr-breadth", type=float, default=0.0061)
    ap.add_argument("--thr-blast", type=float, default=87.0)
    ap.add_argument("--k-required", type=int, default=2)

    ap.add_argument("--bootstrap-n", type=int, default=5000)
    ap.add_argument("--permutation-n", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=450)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metrics, sep="\t")

    # Robust column detection (matches your pipeline naming patterns)
    col_y = pick_col(df, ["y_true", "truth", "label", "truth_cmv_present"])
    col_rpm = pick_col(df, ["log10_rpm_hq", "log10_rpm", "cmv_log10_rpm_mapq", "cmv_log10_rpm"])
    col_breadth = pick_col(df, ["cmv_breadth_mapq", "cmv_breadth", "cmv_breadth_from_cov"])
    col_blast = pick_col(df, ["blast_cmv_mean_identity_top5", "blast_mean_identity_top5", "blast_mean_identity"])

    use = df[[col_y, col_rpm, col_breadth, col_blast]].copy()
    use[col_y] = use[col_y].astype(int)
    for c in [col_rpm, col_breadth, col_blast]:
        use[c] = pd.to_numeric(use[c], errors="coerce")

    # Load thresholds from JSON if provided, otherwise use CLI/defaults
    thr_defaults = (args.thr_log10_rpm, args.thr_breadth, args.thr_blast, args.k_required)
    if args.thresholds_json:
        thr_rpm, thr_breadth, thr_blast, k_required = load_thresholds_json(args.thresholds_json, thr_defaults)
    else:
        thr_rpm, thr_breadth, thr_blast, k_required = thr_defaults

    # Observed performance
    y_true = use[col_y].values
    y_pred = apply_rule(use, col_rpm, col_breadth, col_blast, thr_rpm, thr_breadth, thr_blast, k_required).values
    tp, tn, fp, fn = confusion_counts(y_true, y_pred)
    obs_recall, obs_spec, obs_prec, obs_f1, obs_bal = metrics_from_counts(tp, tn, fp, fn)

    # --------------------------
    # Bootstrap: stability of metrics under resampling
    # --------------------------
    n = len(use)
    boot_recall = np.empty(args.bootstrap_n)
    boot_spec = np.empty(args.bootstrap_n)
    boot_prec = np.empty(args.bootstrap_n)
    boot_f1 = np.empty(args.bootstrap_n)
    boot_bal = np.empty(args.bootstrap_n)

    idx_all = np.arange(n)

    for i in range(args.bootstrap_n):
        idx = rng.choice(idx_all, size=n, replace=True)
        samp = use.iloc[idx]
        yt = samp[col_y].values
        yp = apply_rule(samp, col_rpm, col_breadth, col_blast, thr_rpm, thr_breadth, thr_blast, k_required).values
        tp_b, tn_b, fp_b, fn_b = confusion_counts(yt, yp)
        r, s, p, f, b = metrics_from_counts(tp_b, tn_b, fp_b, fn_b)
        boot_recall[i] = r
        boot_spec[i] = s
        boot_prec[i] = p
        boot_f1[i] = f
        boot_bal[i] = b

    # --------------------------
    # Permutation: significance vs chance
    # --------------------------
    perm_f1 = np.empty(args.permutation_n)
    for i in range(args.permutation_n):
        y_perm = rng.permutation(y_true)
        tp_p, tn_p, fp_p, fn_p = confusion_counts(y_perm, y_pred)  # fixed predictions, shuffled truth
        _, _, _, f, _ = metrics_from_counts(tp_p, tn_p, fp_p, fn_p)
        perm_f1[i] = f

    # one-sided p-value: P(F1_perm >= F1_obs)
    p_value = (np.sum(perm_f1 >= obs_f1) + 1) / (len(perm_f1) + 1)

    # --------------------------
    # Write summary TSV (nice for Results)
    # --------------------------
    def fmt_ci(arr):
        lo, hi = percentile_ci(arr, 2.5, 97.5)
        return lo, hi

    r_lo, r_hi = fmt_ci(boot_recall)
    s_lo, s_hi = fmt_ci(boot_spec)
    p_lo, p_hi = fmt_ci(boot_prec)
    f_lo, f_hi = fmt_ci(boot_f1)
    b_lo, b_hi = fmt_ci(boot_bal)

    summary = pd.DataFrame([{
        "thr_log10_rpm_hq": thr_rpm,
        "thr_breadth": thr_breadth,
        "thr_blast_identity": thr_blast,
        "k_required": k_required,

        "TP": tp, "TN": tn, "FP": fp, "FN": fn,

        "recall_obs": obs_recall, "recall_boot_ci95_lo": r_lo, "recall_boot_ci95_hi": r_hi,
        "specificity_obs": obs_spec, "specificity_boot_ci95_lo": s_lo, "specificity_boot_ci95_hi": s_hi,
        "precision_obs": obs_prec, "precision_boot_ci95_lo": p_lo, "precision_boot_ci95_hi": p_hi,
        "f1_obs": obs_f1, "f1_boot_ci95_lo": f_lo, "f1_boot_ci95_hi": f_hi,
        "balanced_acc_obs": obs_bal, "balanced_acc_boot_ci95_lo": b_lo, "balanced_acc_boot_ci95_hi": b_hi,

        "permutation_n": args.permutation_n,
        "bootstrap_n": args.bootstrap_n,
        "perm_test_stat": "F1",
        "perm_p_value_one_sided": p_value,
    }])

    summary_path = outdir / "resampling_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    # --------------------------
    # FIG 5: Bootstrap distributions
    # --------------------------
    fig = plt.figure(figsize=(11.5, 7.5))
    panels = [
        ("Recall (Sensitivity)", boot_recall, obs_recall),
        ("Specificity",          boot_spec,   obs_spec),
        ("Precision (PPV)",      boot_prec,   obs_prec),
        ("F1 score",             boot_f1,     obs_f1),
    ]

    for i, (title, arr, obs) in enumerate(panels, start=1):
        ax = fig.add_subplot(2, 2, i)
        ax.hist(arr, bins=40, alpha=0.8)
        ax.axvline(obs, linestyle="--", linewidth=2)
        lo, hi = percentile_ci(arr, 2.5, 97.5)
        ax.axvline(lo, linestyle=":", linewidth=2)
        ax.axvline(hi, linestyle=":", linewidth=2)
        ax.set_title(title)
        ax.set_xlabel("Value")
        ax.set_ylabel("Bootstrap count")
        ax.text(
            0.02, 0.98,
            f"Observed = {obs:.3f}\n95% CI = [{lo:.3f}, {hi:.3f}]",
            transform=ax.transAxes,
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", alpha=0.9),
        )

    fig.suptitle(
        "Bootstrap validation of final CMV detection rule (fixed thresholds)\n"
        f"Rule: ≥{k_required}/3 metrics pass | log10RPM≥{thr_rpm}, breadth≥{thr_breadth}, BLAST≥{thr_blast}%",
        y=0.98
    )
    plt.tight_layout()
    fig5_path = outdir / "Fig5_bootstrap_metric_distributions.png"
    plt.savefig(fig5_path, dpi=args.dpi)
    plt.close(fig)

    # --------------------------
    # FIG 6: Permutation null distribution (F1)
    # --------------------------
    plt.figure(figsize=(7.5, 5.6))
    ax = plt.gca()
    ax.hist(perm_f1, bins=45, alpha=0.85)
    ax.axvline(obs_f1, linestyle="--", linewidth=2)
    ax.set_title("Permutation test: null distribution of F1 (labels shuffled)")
    ax.set_xlabel("F1 under permuted labels")
    ax.set_ylabel("Count")
    ax.text(
        0.02, 0.98,
        f"Observed F1 = {obs_f1:.3f}\nPermutation p (one-sided) = {p_value:.4g}\nN permutations = {args.permutation_n}",
        transform=ax.transAxes,
        va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", alpha=0.9),
    )
    plt.tight_layout()
    fig6_path = outdir / "Fig6_permutation_null_F1.png"
    plt.savefig(fig6_path, dpi=args.dpi)
    plt.close()

    # Print quick summary
    print("=== Observed performance (final fixed rule) ===")
    print(f"Thresholds: log10RPM≥{thr_rpm}, breadth≥{thr_breadth}, BLAST≥{thr_blast} ; k≥{k_required}/3")
    print(f"TP={tp} TN={tn} FP={fp} FN={fn}")
    print(f"Recall={obs_recall:.3f} Specificity={obs_spec:.3f} Precision={obs_prec:.3f} F1={obs_f1:.3f} BalAcc={obs_bal:.3f}")
    print(f"Bootstrap N={args.bootstrap_n} | Permutation N={args.permutation_n} | p={p_value:.4g}")
    print(f"Saved: {fig5_path}")
    print(f"Saved: {fig6_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()
