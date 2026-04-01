#!/usr/bin/env python3
"""
make_kofn_confusion_matrices.py

Create confusion matrices for 1/3, 2/3, and 3/3 CMV detection rules
using three selected metrics:
  - log10_rpm_hq
  - cmv_breadth_mapq
  - blast_cmv_mean_identity_top5

Expected input TSV columns:
  - truth_cmv_present (or another truth column via --truth-col)
  - log10_rpm_hq
  - cmv_breadth_mapq
  - blast_cmv_mean_identity_top5

Missing metric values are treated as no evidence of CMV and set to 0.

Outputs:
  - kofn_rule_summary.tsv
  - confusion_matrix_1of3.png
  - confusion_matrix_2of3.png
  - confusion_matrix_3of3.png
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def safe_div(a, b):
    return float(a) / float(b) if b else float("nan")


def confusion_counts(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)

    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn


def metrics_from_counts(tp, tn, fp, fn):
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    precision = safe_div(tp, tp + fp)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    return recall, specificity, precision, f1


def plot_confusion_matrix(tn, fp, fn, tp, metrics_text, title, outpath):
    mat = np.array([[tn, fp],
                    [fn, tp]])

    plt.figure(figsize=(6.2, 5.6))
    im = plt.imshow(mat, aspect="auto")
    plt.colorbar(im, label="Number of samples")

    plt.xticks([0, 1], ["Pred CMV−", "Pred CMV+"])
    plt.yticks([0, 1], ["True CMV−", "True CMV+"])
    plt.title(title, fontsize=15)

    for (i, j), value in np.ndenumerate(mat):
        plt.text(j, i, str(int(value)), ha="center", va="center", fontsize=18)

    plt.figtext(
        0.5, 0.03, metrics_text,
        ha="center", va="bottom", fontsize=10
    )

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def evaluate_rule(df, k, rpm_col, breadth_col, blast_col,
                  rpm_threshold, breadth_threshold, blast_threshold):
    d = df.copy()

    d["pass_rpm"] = pd.to_numeric(d[rpm_col], errors="coerce") >= rpm_threshold
    d["pass_breadth"] = pd.to_numeric(d[breadth_col], errors="coerce") >= breadth_threshold
    d["pass_blast"] = pd.to_numeric(d[blast_col], errors="coerce") >= blast_threshold

    d["n_pass"] = d[["pass_rpm", "pass_breadth", "pass_blast"]].sum(axis=1)
    d["predicted"] = (d["n_pass"] >= k).astype(int)

    y_true = d["truth_cmv_present"].astype(int).to_numpy()
    y_pred = d["predicted"].astype(int).to_numpy()

    tp, tn, fp, fn = confusion_counts(y_true, y_pred)
    recall, specificity, precision, f1 = metrics_from_counts(tp, tn, fp, fn)

    return {
        "rule": f"{k}/3",
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "recall": recall,
        "specificity": specificity,
        "precision": precision,
        "F1": f1,
        "df": d
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input metrics TSV")
    parser.add_argument("--outdir", required=True, help="Output directory")

    parser.add_argument("--truth-col", default="truth_cmv_present")
    parser.add_argument("--rpm-col", default="log10_rpm_hq")
    parser.add_argument("--breadth-col", default="cmv_breadth_mapq")
    parser.add_argument("--blast-col", default="blast_cmv_mean_identity_top5")

    parser.add_argument("--rpm-threshold", type=float, default=2.001474878361447)
    parser.add_argument("--breadth-threshold", type=float, default=0.006115104860680813)
    parser.add_argument("--blast-threshold", type=float, default=86.9904375)

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")

    required_cols = [
        args.truth_col,
        args.rpm_col,
        args.breadth_col,
        args.blast_col,
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    df = df.rename(columns={args.truth_col: "truth_cmv_present"}).copy()

    df["truth_cmv_present"] = pd.to_numeric(df["truth_cmv_present"], errors="coerce")
    df = df.dropna(subset=["truth_cmv_present"]).copy()
    df["truth_cmv_present"] = df["truth_cmv_present"].astype(int)

    # Treat missing metric values as no evidence of CMV
    df[args.rpm_col] = pd.to_numeric(df[args.rpm_col], errors="coerce").fillna(0)
    df[args.breadth_col] = pd.to_numeric(df[args.breadth_col], errors="coerce").fillna(0)
    df[args.blast_col] = pd.to_numeric(df[args.blast_col], errors="coerce").fillna(0)

    results = []
    all_rule_outputs = {}

    for k in [1, 2, 3]:
        res = evaluate_rule(
            df=df,
            k=k,
            rpm_col=args.rpm_col,
            breadth_col=args.breadth_col,
            blast_col=args.blast_col,
            rpm_threshold=args.rpm_threshold,
            breadth_threshold=args.breadth_threshold,
            blast_threshold=args.blast_threshold,
        )
        results.append({
            "rule": res["rule"],
            "TP": res["TP"],
            "TN": res["TN"],
            "FP": res["FP"],
            "FN": res["FN"],
            "recall": res["recall"],
            "specificity": res["specificity"],
            "precision": res["precision"],
            "F1": res["F1"],
        })
        all_rule_outputs[k] = res

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(outdir / "kofn_rule_summary.tsv", sep="\t", index=False)

    for k in [1, 2, 3]:
        res = all_rule_outputs[k]

        metrics_text = (
            f"TP={res['TP']}  TN={res['TN']}  FP={res['FP']}  FN={res['FN']}\n"
            f"Recall={res['recall']:.3f}  Specificity={res['specificity']:.3f}  "
            f"Precision={res['precision']:.3f}  F1={res['F1']:.3f}"
        )

        plot_confusion_matrix(
            tn=res["TN"],
            fp=res["FP"],
            fn=res["FN"],
            tp=res["TP"],
            metrics_text=metrics_text,
            title=f"Confusion matrix: {k}/3 rule",
            outpath=outdir / f"confusion_matrix_{k}of3.png"
        )

    print("Done")
    print(f"Saved summary: {outdir / 'kofn_rule_summary.tsv'}")
    print("Saved figures:")
    print(f"  {outdir / 'confusion_matrix_1of3.png'}")
    print(f"  {outdir / 'confusion_matrix_2of3.png'}")
    print(f"  {outdir / 'confusion_matrix_3of3.png'}")


if __name__ == "__main__":
    main()