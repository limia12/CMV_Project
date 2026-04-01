#!/usr/bin/env python3
"""
make_final_results_plots_v2.py

Results-ready figures for the FINAL multi-metric CMV detection rule:
  - Clean 3D scatter plot (3 metrics) with subtle threshold planes + TP/TN/FP/FN categories
  - 2D projections (RPM vs breadth, RPM vs BLAST, breadth vs BLAST) with threshold lines
  - Confusion matrix + key metrics

Rule:
  CMV positive if >= k of 3 metrics pass:
    - log10 RPM (HQ CMV reads) >= thr_log10_rpm
    - CMV breadth (MAPQ-filtered) >= thr_breadth
    - BLAST mean % identity (top 5 hits) >= thr_blast

Usage:
  python3 make_final_results_plots_v2.py \
    --metrics /path/to/metrics_table.tsv \
    --outdir /path/to/results_figures \
    --thr-log10-rpm 2.0 --thr-breadth 0.0061 --thr-blast 87 \
    --k-required 2
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------
# Helpers
# ----------------------------
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


def compute_confusion(y_true, y_pred):
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
    prec = safe_div(tp, tp + fp)            # precision / PPV
    f1 = safe_div(2 * prec * recall, prec + recall)
    return recall, spec, prec, f1


def apply_k_of_3_rule(df, col_rpm, col_breadth, col_blast, thr_rpm, thr_breadth, thr_blast, k_required):
    """
    Conservative handling of missing values:
      NA counts as FAIL for that metric (clinical-friendly).
    """
    pass_sum = (
        (df[col_rpm] >= thr_rpm).fillna(False).astype(int)
        + (df[col_breadth] >= thr_breadth).fillna(False).astype(int)
        + (df[col_blast] >= thr_blast).fillna(False).astype(int)
    )
    return (pass_sum >= k_required).astype(int), pass_sum


# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True, help="metrics_table.tsv")
    ap.add_argument("--outdir", required=True, help="Output directory for figures")

    ap.add_argument("--thr-log10-rpm", type=float, default=2.0)
    ap.add_argument("--thr-breadth", type=float, default=0.0061)
    ap.add_argument("--thr-blast", type=float, default=87.0)
    ap.add_argument("--k-required", type=int, default=2)

    ap.add_argument("--dpi", type=int, default=450)
    ap.add_argument("--elev", type=float, default=22.0, help="3D view elevation")
    ap.add_argument("--azim", type=float, default=-62.0, help="3D view azimuth")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metrics, sep="\t")

    # Robust column detection (matches your pipeline)
    col_y = pick_col(df, ["y_true", "truth", "label", "truth_cmv_present"])
    col_rpm = pick_col(df, ["log10_rpm_hq", "log10_rpm", "cmv_log10_rpm_mapq", "cmv_log10_rpm"])
    col_breadth = pick_col(df, ["cmv_breadth_mapq", "cmv_breadth", "cmv_breadth_from_cov"])
    col_blast = pick_col(df, ["blast_cmv_mean_identity_top5", "blast_mean_identity_top5", "blast_mean_identity"])

    use = df[[col_y, col_rpm, col_breadth, col_blast]].copy()
    use[col_y] = use[col_y].astype(int)
    for c in [col_rpm, col_breadth, col_blast]:
        use[c] = pd.to_numeric(use[c], errors="coerce")

    # Apply final rule
    y_pred, pass_sum = apply_k_of_3_rule(
        use, col_rpm, col_breadth, col_blast,
        args.thr_log10_rpm, args.thr_breadth, args.thr_blast, args.k_required
    )
    y_true = use[col_y].values

    # Confusion + metrics
    tp, tn, fp, fn = compute_confusion(y_true, y_pred.values)
    recall, spec, prec, f1 = metrics_from_counts(tp, tn, fp, fn)

    # Category per sample
    cat = np.full(len(use), "", dtype=object)
    cat[(y_true == 1) & (y_pred.values == 1)] = "TP"
    cat[(y_true == 0) & (y_pred.values == 0)] = "TN"
    cat[(y_true == 0) & (y_pred.values == 1)] = "FP"
    cat[(y_true == 1) & (y_pred.values == 0)] = "FN"

    # Add pass count for optional point sizing
    use["_pass_sum"] = pass_sum.values
    use["_cat"] = cat

    # ----------------------------
    # Figure 1: Cleaner 3D scatter
    # ----------------------------
    fig = plt.figure(figsize=(9.2, 7.2))
    ax = fig.add_subplot(111, projection="3d")

    # Reduce 3D clutter
    ax.grid(False)
    try:
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
    except Exception:
        pass

    # Use BOTH marker + colour to reduce confusion
    # (Colours are standard Matplotlib “tab:” palette; clear for readers)
    style = {
        "TP": dict(marker="o", color="tab:blue"),
        "TN": dict(marker="s", color="tab:orange"),
        "FP": dict(marker="x", color="tab:green"),
        "FN": dict(marker="^", color="tab:red"),
    }

    # Point size can encode “how many metrics passed” (nice explanatory touch)
    # scale: 1..3 -> 35..75
    sizes = 35 + 20 * use["_pass_sum"].clip(0, 3)

    for k in ["TP", "TN", "FP", "FN"]:
        m = (use["_cat"] == k)
        if not np.any(m):
            continue
        ax.scatter(
            use.loc[m, col_rpm],
            use.loc[m, col_breadth],
            use.loc[m, col_blast],
            s=sizes[m],
            alpha=0.85,
            linewidths=0.7,
            label=f"{k} (n={int(m.sum())})",
            **style[k],
        )

    # Subtle threshold planes (filled surfaces, low alpha) instead of wireframes
    x_min, x_max = float(np.nanmin(use[col_rpm])), float(np.nanmax(use[col_rpm]))
    y_min, y_max = float(np.nanmin(use[col_breadth])), float(np.nanmax(use[col_breadth]))
    z_min, z_max = float(np.nanmin(use[col_blast])), float(np.nanmax(use[col_blast]))

    # plane x = thr_rpm
    Y, Z = np.meshgrid(np.linspace(y_min, y_max, 2), np.linspace(z_min, z_max, 2))
    X = np.full_like(Y, args.thr_log10_rpm)
    ax.plot_surface(X, Y, Z, alpha=0.08, linewidth=0)

    # plane y = thr_breadth
    X2, Z2 = np.meshgrid(np.linspace(x_min, x_max, 2), np.linspace(z_min, z_max, 2))
    Y2 = np.full_like(X2, args.thr_breadth)
    ax.plot_surface(X2, Y2, Z2, alpha=0.08, linewidth=0)

    # plane z = thr_blast
    X3, Y3 = np.meshgrid(np.linspace(x_min, x_max, 2), np.linspace(y_min, y_max, 2))
    Z3 = np.full_like(X3, args.thr_blast)
    ax.plot_surface(X3, Y3, Z3, alpha=0.08, linewidth=0)

    ax.set_xlabel("log10 CMV RPM (HQ)")
    ax.set_ylabel("CMV breadth (MAPQ-filtered)")
    ax.set_zlabel("BLAST CMV % identity (top 5 mean)")

    ax.view_init(elev=args.elev, azim=args.azim)

    ax.set_title(
        f"3-metric CMV detection space (≥{args.k_required}/3 rule)\n"
        f"Thresholds: log10RPM≥{args.thr_log10_rpm}, breadth≥{args.thr_breadth}, BLAST≥{args.thr_blast}%\n"
        f"Point size = #metrics passing (1–3)"
    )

    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), frameon=True)
    plt.tight_layout()
    fig1_path = outdir / "Fig1_3D_metrics_TP_TN_FP_FN_clean.png"
    plt.savefig(fig1_path, dpi=args.dpi)
    plt.close(fig)

    # ----------------------------
    # Figure 2: 2D projection panel (much easier for readers)
    # ----------------------------
    fig = plt.figure(figsize=(12.0, 4.2))

    pairs = [
        (col_rpm, col_breadth, "log10 CMV RPM (HQ)", "CMV breadth (MAPQ-filtered)", args.thr_log10_rpm, args.thr_breadth),
        (col_rpm, col_blast, "log10 CMV RPM (HQ)", "BLAST CMV % identity (top 5 mean)", args.thr_log10_rpm, args.thr_blast),
        (col_breadth, col_blast, "CMV breadth (MAPQ-filtered)", "BLAST CMV % identity (top 5 mean)", args.thr_breadth, args.thr_blast),
    ]

    for i, (cx, cy, xl, yl, tx, ty) in enumerate(pairs, start=1):
        ax = fig.add_subplot(1, 3, i)
        for k in ["TP", "TN", "FP", "FN"]:
            m = (use["_cat"] == k)
            if not np.any(m):
                continue
            ax.scatter(
                use.loc[m, cx],
                use.loc[m, cy],
                s=35,
                alpha=0.85,
                linewidths=0.7,
                label=k if i == 1 else None,
                **style[k],
            )
        ax.axvline(tx, linestyle="--", linewidth=1.5)
        ax.axhline(ty, linestyle="--", linewidth=1.5)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title("")

    fig.suptitle(
        f"2D projections of the 3-metric CMV rule (≥{args.k_required}/3)\n"
        f"Dashed lines = thresholds",
        y=1.04
    )
    # single legend
    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=True)
    plt.tight_layout()
    fig2_path = outdir / "Fig2_2D_projections_TP_TN_FP_FN.png"
    plt.savefig(fig2_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    # ----------------------------
    # Figure 3: Confusion matrix
    # ----------------------------
    cm = np.array([[tn, fp],
                   [fn, tp]])

    fig = plt.figure(figsize=(6.3, 5.7))
    ax = fig.add_subplot(111)
    ax.imshow(cm)

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Pred CMV−", "Pred CMV+"])
    ax.set_yticklabels(["True CMV−", "True CMV+"])
    ax.set_title("Confusion matrix (final multi-metric rule)")

    for (i, j), v in np.ndenumerate(cm):
        ax.text(j, i, str(int(v)), ha="center", va="center", fontsize=16)

    summary = (
        f"TP={tp}  TN={tn}  FP={fp}  FN={fn}\n"
        f"Recall={recall:.3f}  Specificity={spec:.3f}  Precision={prec:.3f}  F1={f1:.3f}"
    )
    ax.set_xlabel(summary)

    plt.tight_layout()
    fig3_path = outdir / "Fig3_confusion_matrix.png"
    plt.savefig(fig3_path, dpi=args.dpi)
    plt.close(fig)

    # Print summary
    print("=== Final rule summary (computed from metrics_table.tsv) ===")
    print(f"Thresholds: log10RPM≥{args.thr_log10_rpm}, breadth≥{args.thr_breadth}, BLAST≥{args.thr_blast} ; k≥{args.k_required}/3")
    print(f"TP={tp} TN={tn} FP={fp} FN={fn}")
    print(f"Recall={recall:.3f} Specificity={spec:.3f} Precision={prec:.3f} F1={f1:.3f}")
    print(f"Saved:\n  {fig1_path}\n  {fig2_path}\n  {fig3_path}")


if __name__ == "__main__":
    main()
