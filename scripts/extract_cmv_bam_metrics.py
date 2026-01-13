#!/usr/bin/env python3
"""
bam_metrics.py

This module calculates simple coverage metrics from CMV-aligned BAM files.

It produces three CSV outputs
per sample:

1) <sample>_summary.csv
   - One row per sample
   - Includes total read count and a simple "average genome depth" value

2) <sample>_gene.csv
   - One row per gene region per sample
   - Includes average depth and breadth of coverage for each gene region

3) <sample>_per_base.csv
   - One row per base position per gene region per sample
   - Includes per-base depth values (coverage at each position)

Important notes 
----------------
- This script uses hardcoded CMV Merlin genome gene coordinates.

"""

import os

import pandas as pd
import pysam

# Hardcoded gene regions for CMV Merlin genome
#
# Each tuple is:
#   (contig_name, start_position, end_position, gene_name)
#
# These coordinates are used exactly as provided. In pysam, start is
# 0-based and end is exclusive for count_coverage. That means the region
# covers positions start .. end-1.
CMV_GENES = [
    ("merlinGenome", 172328, 174090, "IE1"),
    ("merlinGenome", 78193, 81922, "UL54"),
    ("merlinGenome", 82065, 84789, "gpUL55"),
    ("merlinGenome", 109223, 111452, "UL75"),
    ("merlinGenome", 141797, 143921, "UL97"),
    ("merlinGenome", 120655, 122341, "UL83"),
]


def analyze_bam(bam_path):
    """
    Analyze one BAM file and calculate summary, per-gene, and per-base metrics.

    Parameters
    ----------
    bam_path : str
        Path to a BAM file aligned to the CMV Merlin reference.

    Returns
    -------
    tuple
        (sample, summary, gene_stats, per_base_depth_records)

        sample : str
            Sample name derived from the BAM filename (basename without .bam).

        summary : dict
            A one-row summary for the sample:
              - sample: sample name
              - total_reads: total reads counted from the BAM
              - avg_genome_depth: average of per-gene average depths

            Note:
            "avg_genome_depth" is not a true whole-genome metric. It is the
            average of the average depths across the hardcoded gene regions.

        gene_stats : list of dict
            A list of per-gene results. Each dict contains:
              - sample: sample name
              - gene: gene name
              - avg_depth: mean depth across that gene region
              - breadth: fraction of bases with depth > 0 across that region

        per_base_depth_records : list of dict
            A list of per-base depth values. Each dict contains:
              - sample: sample name
              - gene: gene name
              - position: genome coordinate for that base
              - depth: depth at that base

    Step-by-step explanation (what happens inside this function)
    ------------------------------------------------------------
    1) Derive a sample name from the BAM filename.
    2) Open the BAM file using pysam.
    3) Count total reads in the BAM.
    4) For each gene region in CMV_GENES:
       a) Use bam.count_coverage() to get base-by-base coverage.
       b) Convert A/C/G/T coverage into a single depth per position by summing.
       c) Compute average depth and breadth of coverage for the region.
       d) Store per-base depth records for later CSV output.
       e) If anything fails (e.g., contig name mismatch), record zeros.
    5) Compute avg_genome_depth as the mean of all per-gene avg_depth values.
    6) Return everything needed for downstream CSV writing.
    """
    sample = os.path.basename(bam_path).replace(".bam", "")
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Total reads in BAM. This is a quick overall size/volume metric.
    total_reads = bam.count()

    # gene_stats will hold one row per gene region for this sample.
    gene_stats = []

    # per_base_depth_records will hold one row per base per gene region.
    per_base_depth_records = []

    for chrom, start, end, gene in CMV_GENES:
        try:
            # count_coverage returns 4 arrays: A, C, G, T counts per position.
            # We sum them to get total depth per base.
            coverage = bam.count_coverage(contig=chrom, start=start, end=end)
            depths = [sum(pos) for pos in zip(*coverage)]

            # Region length. Used for averaging and for breadth calculations.
            length = end - start

            # Total depth is the sum of per-base depths across the region.
            total_depth = sum(depths)

            # Average depth: total_depth divided by number of bases in region.
            avg_depth = total_depth / length if length else 0

            # Breadth: fraction of bases with depth > 0.
            covered_bases = sum(1 for d in depths if d > 0)
            breadth = covered_bases / length if length else 0

            # Store per-base depths with genome coordinates.
            for i, depth in enumerate(depths):
                genome_pos = start + i
                per_base_depth_records.append(
                    {
                        "sample": sample,
                        "gene": gene,
                        "position": genome_pos,
                        "depth": depth,
                    }
                )

        except Exception:
            # If the region cannot be processed (e.g., contig name not found),
            # record zeros so the output stays consistent.
            avg_depth = 0
            breadth = 0

        # Store per-gene summary values.
        gene_stats.append(
            {
                "sample": sample,
                "gene": gene,
                "avg_depth": round(avg_depth, 2),
                "breadth": round(breadth, 3),
            }
        )

    # A simple "average genome depth" based on the mean of the per-gene
    # average depths. This is useful as a single-number summary, but it is
    # not a true whole-genome depth value.
    avg_genome_depth = (
        round(sum(g["avg_depth"] for g in gene_stats) / len(gene_stats), 2)
        if gene_stats
        else 0
    )

    summary = {
        "sample": sample,
        "total_reads": total_reads,
        "avg_genome_depth": avg_genome_depth,
    }

    return sample, summary, gene_stats, per_base_depth_records


def run_bam_metrics(bam_dir, out_dir, log_func=print):
    """
    Run coverage metrics for every BAM file in a directory and write CSVs.

    Parameters
    ----------
    bam_dir : str
        Directory containing one or more BAM files. This function scans for
        files ending in ".bam".

    out_dir : str
        Directory where output CSV files will be written.

    log_func : callable, optional
        A logging function used to print progress updates. Defaults to print.
        This is helpful if you want to plug this into a GUI (e.g., Streamlit)
        and replace print() with a UI logger.

    Returns
    -------
    bool
        True if processing ran (and at least one BAM existed), otherwise False.

    Step-by-step explanation (what happens inside this function)
    ------------------------------------------------------------
    1) Create the output directory if it does not exist.
    2) Find all BAM files in bam_dir.
    3) If there are no BAM files, log an error message and return False.
    4) For each BAM file:
       a) Call analyze_bam() to compute metrics.
       b) Create three output paths for summary, per-gene, and per-base results.
       c) Write the results to CSV files using pandas.
       d) Log where outputs were saved.
    5) Return True when complete.
    """
    os.makedirs(out_dir, exist_ok=True)
    bam_files = sorted(f for f in os.listdir(bam_dir) if f.endswith(".bam"))

    if not bam_files:
        log_func("❌ No BAM files found in the specified directory.")
        return False

    for bam_file in bam_files:
        full_path = os.path.join(bam_dir, bam_file)
        log_func(f"🔍 Processing {bam_file}")
        sample, summary, gene_data, per_base_data = analyze_bam(full_path)

        # Output files for this sample.
        summary_path = os.path.join(out_dir, f"{sample}_summary.csv")
        gene_path = os.path.join(out_dir, f"{sample}_gene.csv")
        per_base_path = os.path.join(out_dir, f"{sample}_per_base.csv")

        # Write CSV outputs.
        pd.DataFrame([summary]).to_csv(summary_path, index=False)
        pd.DataFrame(gene_data).to_csv(gene_path, index=False)
        pd.DataFrame(per_base_data).to_csv(per_base_path, index=False)

        log_func(
            "✅ Outputs saved for sample "
            f"'{sample}' to:\n  - {summary_path}\n  - {gene_path}\n  - "
            f"{per_base_path}"
        )

    return True
