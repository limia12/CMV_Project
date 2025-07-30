#!/usr/bin/env python3

import os
import pysam
import pandas as pd

# Hardcoded gene regions for CMV Merlin genome
CMV_GENES = [
    ("merlinGenome", 172328, 174090, "IE1"),
    ("merlinGenome", 78193, 81922, "UL54"),
    ("merlinGenome", 82065, 84789, "gpUL55"),
    ("merlinGenome", 109223, 111452, "UL75"),
    ("merlinGenome", 141797, 143921, "UL97"),
    ("merlinGenome", 120655, 122341, "UL83"),
]

def analyze_bam(bam_path):
    sample = os.path.basename(bam_path).replace(".bam", "")
    bam = pysam.AlignmentFile(bam_path, "rb")
    total_reads = bam.count()
    gene_stats = []
    per_base_depth_records = []

    for chrom, start, end, gene in CMV_GENES:
        try:
            coverage = bam.count_coverage(contig=chrom, start=start, end=end)
            depths = [sum(pos) for pos in zip(*coverage)]
            length = end - start
            total_depth = sum(depths)
            avg_depth = total_depth / length if length else 0
            covered_bases = sum(1 for d in depths if d > 0)
            breadth = covered_bases / length if length else 0

            for i, depth in enumerate(depths):
                genome_pos = start + i
                per_base_depth_records.append({
                    "sample": sample,
                    "gene": gene,
                    "position": genome_pos,
                    "depth": depth
                })

        except Exception:
            avg_depth = 0
            breadth = 0

        gene_stats.append({
            "sample": sample,
            "gene": gene,
            "avg_depth": round(avg_depth, 2),
            "breadth": round(breadth, 3)
        })

    avg_genome_depth = round(
        sum(g["avg_depth"] for g in gene_stats) / len(gene_stats), 2
    ) if gene_stats else 0

    summary = {
        "sample": sample,
        "total_reads": total_reads,
        "avg_genome_depth": avg_genome_depth
    }

    return sample, summary, gene_stats, per_base_depth_records

def run_bam_metrics(bam_dir, out_dir, log_func=print):
    os.makedirs(out_dir, exist_ok=True)
    bam_files = sorted(f for f in os.listdir(bam_dir) if f.endswith(".bam"))

    if not bam_files:
        log_func("❌ No BAM files found in the specified directory.")
        return False

    for bam_file in bam_files:
        full_path = os.path.join(bam_dir, bam_file)
        log_func(f"🔍 Processing {bam_file}")
        sample, summary, gene_data, per_base_data = analyze_bam(full_path)

        # Output files
        summary_path = os.path.join(out_dir, f"{sample}_summary.csv")
        gene_path = os.path.join(out_dir, f"{sample}_gene.csv")
        per_base_path = os.path.join(out_dir, f"{sample}_per_base.csv")

        pd.DataFrame([summary]).to_csv(summary_path, index=False)
        pd.DataFrame(gene_data).to_csv(gene_path, index=False)
        pd.DataFrame(per_base_data).to_csv(per_base_path, index=False)

        log_func(f"✅ Outputs saved for sample '{sample}' to:\n  - {summary_path}\n  - {gene_path}\n  - {per_base_path}")
    
    return True
