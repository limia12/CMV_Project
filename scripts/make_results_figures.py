#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of {candidates} found. Available: {list(df.columns)[:80]}")

def confusion(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn

def safe_div(a, b):
    return a / b if b else 0.0

def metrics_from_cm(tp, tn, fp, fn):
    recall = safe_div(tp, tp + fn)
    spec = safe_div(tn, tn + fp)
    prec = safe_div(tp, tp + fp)
    f1 = safe_div(2 * prec * recall, prec + recall)
    return recall, spec, prec, f1

def sweep_threshold(metric_values, y_true):
    v = np.asarray(metric_values, dtype=float)
    y = np.asarray(y_true).astype(int)
    mask = ~np.isnan(v)
    v = v[mask]; y = y[mask]
    thr_list = np.unique(v)
    thr_list.sort()

    recalls, specs, precs, f1s = [], [], [], []
    for thr in thr_list:
        yp = (v >= thr).astype(int)
        tp, tn, fp, fn = confusion(y, yp)
        r, s, p, f = metrics_from_cm(tp, tn, fp, fn)
        recalls.append(r); specs.append(s); precs.append(p); f1s.append(f)
    return thr_list, np.array(recalls), np.array(specs), np.array(precs), np.array(f1s)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True, help="metrics_table.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--thresholds-json", default="", help="optional thresholds.json")
    ap.add_argument("--thr-log10-rpm", type=float, default=2.0)
    ap.add_argument("--thr-breadth", type=float, default=0.0061)
    ap.add_argument("--thr-blast", type=float, default=87.0)
    ap.add_argument("--k-required", type=int, default=2)
    ap.add_argument("--dpi", type=int, default=450)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metrics, sep="\t")

    col_y = pick_col(df, ["y_true", "truth", "label", "truth_cmv_present"])
    col_rpm = pick_col(df, ["log10_rpm_hq", "log10_rpm", "cmv_log10_rpm_mapq", "cmv_log10_rpm"])
    col_breadth = pick_col(df, ["cmv_breadth_mapq", "cmv_breadth", "cmv_breadth_from_cov"])
    col_blast = pick_col(df, ["blast_cmv_mean_identity_top5", "blast_mean_identity_top5", "blast_mean_identity"])

    use = df[[col_y, col_rpm, col_breadth, col_blast]].copy()
    use[col_y] = use[col_y].astype(int)
    for c in [col_rpm, col_breadth, col_blast]:
        use[c] = pd.to_numeric(use[c], errors="coerce")

    thr_rpm, thr_breadth, thr_blast, k_required = args.thr_log10_rpm, args.thr_breadth, args.thr_blast, args.k_required

    if args.thresholds_json:
        obj = json.loads(Path(args.thresholds_json).read_text())
        if isinstance(obj, dict):
            thr_rpm = float(obj.get("thr_log10_rpm_hq", obj.get("thr_rpm", thr_rpm)))
            thr_breadth = float(obj.get("thr_cmv_breadth_mapq", obj.get("thr_breadth", thr_breadth)))
            thr_blast = float(obj.get("thr_blast_identity", obj.get("thr_blast", thr_blast)))
            k_required = int(obj.get("k_required", obj.get("k", k_required)))

    y_true = use[col_y].values

    pass_sum = (
        (use[col_rpm] >= thr_rpm).fillna(False).astype(int)
        + (use[col_breadth] >= thr_breadth).fillna(False).astype(int)
        + (use[col_blast] >= thr_blast).fillna(False).astype(int)
    )
    y_pred = (pass_sum >= k_required).astype(int)

    tp, tn, fp, fn = confusion(y_true, y_pred)
    recall, spec, prec, f1 = metrics_from_cm(tp, tn, fp, fn)

    # Categories TP/TN/FP/FN for later plots
    cat = np.full(len(use), "", dtype=object)
    cat[(y_true == 1) & (y_pred == 1)] = "TP"
    cat[(y_true == 0) & (y_pred == 0)] = "TN"
    cat[(y_true == 0) & (y_pred == 1)] = "FP"
    cat[(y_true == 1) & (y_pred == 0)] = "FN"

    # ---- Fig1: confusion matrix ----
    cm = np.array([[tn, fp],
                   [fn, tp]])

    plt.figure(figsize=(6.0, 5.4))
    ax = plt.gca()
    ax.imshow(cm)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["Pred CMV−", "Pred CMV+"])
    ax.set_yticklabels(["True CMV−", "True CMV+"])
    ax.set_title("Confusion matrix (final multi-metric threshold)")
    for (i, j), v in np.ndenumerate(cm):
        ax.text(j, i, str(int(v)), ha="center", va="center", fontsize=16)
    ax.set_xlabel(
        f"TP={tp}  TN={tn}  FP={fp}  FN={fn}\n"
        f"Recall={recall:.3f}  Specificity={spec:.3f}  Precision={prec:.3f}  F1={f1:.3f}"
    )
    plt.tight_layout()
    plt.savefig(outdir / "Fig1_confusion_matrix.png", dpi=args.dpi)
    plt.close()

    # ---- Fig2: distributions + threshold lines ----
    fig = plt.figure(figsize=(10.5, 6.5))
    metrics = [
        (col_rpm, "log10 RPM (HQ CMV reads)", thr_rpm),
        (col_breadth, "CMV genome breadth (MAPQ-filtered)", thr_breadth),
        (col_blast, "BLAST CMV % identity (top 5 mean)", thr_blast),
    ]
    for i, (col, label, thr) in enumerate(metrics, start=1):
        ax = fig.add_subplot(2, 2, i)
        pos = use[use[col_y] == 1][col].dropna().values
        neg = use[use[col_y] == 0][col].dropna().values
        ax.hist(neg, bins=30, alpha=0.5, label="True CMV−")
        ax.hist(pos, bins=30, alpha=0.5, label="True CMV+")
        ax.axvline(thr, linestyle="--", linewidth=2)
        ax.set_title(label)
        ax.set_ylabel("Count")
        ax.legend(frameon=False)
    fig.add_subplot(2, 2, 4).axis("off")
    fig.suptitle("Distributions of candidate metrics with selected thresholds", y=0.98)
    plt.tight_layout()
    plt.savefig(outdir / "Fig2_metric_distributions.png", dpi=args.dpi)
    plt.close(fig)

    # ---- Fig3: trade-off curves ----
    fig = plt.figure(figsize=(12, 8))
    for i, (col, label, thr_chosen) in enumerate(metrics, start=1):
        ax = fig.add_subplot(2, 2, i)
        thr_list, r, s, p, f = sweep_threshold(use[col].values, y_true)
        ax.plot(thr_list, f, linewidth=2, label="F1")
        ax.plot(thr_list, r, linewidth=1.5, label="Recall")
        ax.plot(thr_list, s, linewidth=1.5, label="Specificity")
        ax.axvline(thr_chosen, linestyle="--", linewidth=2)
        if len(f):
            best_idx = int(np.nanargmax(f))
            ax.scatter([thr_list[best_idx]], [f[best_idx]], s=60)
            ax.text(thr_list[best_idx], f[best_idx], "  best F1", va="center")
        ax.set_title(f"{label}\n(single-metric classifier)")
        ax.set_xlabel("Threshold")
        ax.set_ylabel("Score")
        ax.set_ylim(-0.02, 1.02)
        ax.legend(frameon=False)
    fig.add_subplot(2, 2, 4).axis("off")
    fig.suptitle("Trade-off curves justify the selected threshold values", y=0.98)
    plt.tight_layout()
    plt.savefig(outdir / "Fig3_threshold_tradeoffs.png", dpi=args.dpi)
    plt.close(fig)

    # ---- Fig4: decision regions (RPM vs breadth) with BLAST pass/fail ----
    x_min, x_max = float(use[col_rpm].min()), float(use[col_rpm].max())
    y_min, y_max = float(use[col_breadth].min()), float(use[col_breadth].max())
    xs = np.linspace(x_min, x_max, 200)
    ys = np.linspace(y_min, y_max, 200)
    X, Y = np.meshgrid(xs, ys)

    def predict_grid(blast_pass):
        rpm_pass = (X >= thr_rpm).astype(int)
        br_pass = (Y >= thr_breadth).astype(int)
        bl_pass = 1 if blast_pass else 0
        ps = rpm_pass + br_pass + bl_pass
        return (ps >= k_required).astype(int)

    Z_lo = predict_grid(False)
    Z_hi = predict_grid(True)

    fig = plt.figure(figsize=(12, 5.5))
    for idx, (Z, title) in enumerate([
        (Z_lo, f"BLAST identity < {thr_blast}% (fails BLAST)"),
        (Z_hi, f"BLAST identity ≥ {thr_blast}% (passes BLAST)")
    ], start=1):
        ax = fig.add_subplot(1, 2, idx)
        ax.imshow(Z, origin="lower", extent=[x_min, x_max, y_min, y_max], aspect="auto", alpha=0.25)
        for lab, marker in [("TP", "o"), ("TN", "s"), ("FP", "x"), ("FN", "^")]:
            m = (cat == lab)
            if np.any(m):
                ax.scatter(use.loc[m, col_rpm], use.loc[m, col_breadth], marker=marker, s=35, alpha=0.85, label=lab)
        ax.axvline(thr_rpm, linestyle="--", linewidth=2)
        ax.axhline(thr_breadth, linestyle="--", linewidth=2)
        ax.set_xlabel("log10 RPM (HQ)")
        ax.set_ylabel("CMV breadth")
        ax.set_title(title)
        ax.legend(frameon=False)
    fig.suptitle("Why a ≥2-of-3 rule works: decision regions in (RPM, breadth) space", y=0.98)
    plt.tight_layout()
    plt.savefig(outdir / "Fig4_decision_regions.png", dpi=args.dpi)
    plt.close(fig)

    print("=== Final rule summary ===")
    print(f"Thresholds: log10RPM>={thr_rpm}, breadth>={thr_breadth}, BLAST>={thr_blast}; k>={k_required}/3")
    print(f"TP={tp} TN={tn} FP={fp} FN={fn}")
    print(f"Recall={recall:.3f} Specificity={spec:.3f} Precision={prec:.3f} F1={f1:.3f}")
    print(f"Saved figures to: {outdir}")


if __name__ == "__main__":
    main()
