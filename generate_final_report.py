import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- Conserved CMV Regions ---
CONSERVE_REGIONS = [
    {"name": "IE1", "start": 172328, "end": 174090},
    {"name": "UL54", "start": 78193, "end": 81922},
    {"name": "gpUL55", "start": 82065, "end": 84789},
    {"name": "UL75", "start": 109223, "end": 111452},
    {"name": "UL97", "start": 141797, "end": 143921},
    {"name": "UL83", "start": 120655, "end": 122341},
]

# ---------- NEW: helpers to auto-load & merge all samples ----------

def _read_csv_list(files, **read_kwargs):
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, **read_kwargs)
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            print(f"❌ Error reading {f}: {e}")
    if dfs:
        out = pd.concat(dfs, ignore_index=True)
        # standardize dtypes where appropriate
        for col in ("avg_depth", "breadth", "depth", "position", "total_reads"):
            if col in out.columns:
                out[col] = pd.to_numeric(out[col], errors="coerce")
        return out
    return pd.DataFrame()

def _load_all_metrics_from_disk(out_dir, blast_dir, passed_gene_df, passed_blast_df, passed_bam_df, passed_per_base_df):
    """
    Keep the API unchanged but prefer disk-merged data if present.
    Fall back to the passed dataframes only if no disk files exist.
    """
    # Gene-level coverage
    gene_files = glob.glob(os.path.join(out_dir, "*_gene.csv"))
    gene_df_disk = _read_csv_list(gene_files)
    gene_df = gene_df_disk if not gene_df_disk.empty else (passed_gene_df.copy() if passed_gene_df is not None else pd.DataFrame())

    # Per-base depth
    per_base_files = glob.glob(os.path.join(out_dir, "*_per_base.csv"))
    per_base_df_disk = _read_csv_list(per_base_files)
    per_base_df = per_base_df_disk if not per_base_df_disk.empty else (passed_per_base_df.copy() if passed_per_base_df is not None else pd.DataFrame())

    # BAM summaries
    summary_files = glob.glob(os.path.join(out_dir, "*_summary.csv"))
    bam_df_disk = _read_csv_list(summary_files)
    bam_df = bam_df_disk if not bam_df_disk.empty else (passed_bam_df.copy() if passed_bam_df is not None else pd.DataFrame())

    # BLAST top-5 summaries (already per-sample output by summarize_blast_hits.py)
    blast_top5_files = glob.glob(os.path.join(out_dir, "*_blast_top5.csv"))
    blast_df_disk = _read_csv_list(blast_top5_files)
    blast_df = blast_df_disk if not blast_df_disk.empty else (passed_blast_df.copy() if passed_blast_df is not None else pd.DataFrame())

    # Light sanity prints to help debugging
    def _uniq(col, df):
        return list(df[col].dropna().unique()) if col in df.columns and not df.empty else []

    print(f"ℹ️ gene_df samples:    {_uniq('sample', gene_df)}")
    print(f"ℹ️ per_base_df samples:{_uniq('sample', per_base_df)}")
    print(f"ℹ️ bam_df samples:     {_uniq('sample', bam_df)}")
    print(f"ℹ️ blast_df samples:   {_uniq('sample', blast_df)}")

    return gene_df, blast_df, bam_df, per_base_df

# -------------------------------------------------------------------

# --- BLAST identity heatmap ---
def plot_blast_identity_heatmaps(blast_dir, out_dir, min_identity=85, bin_size=100):
    blast_files = glob.glob(os.path.join(blast_dir, "*_blast.tsv"))
    if not blast_files:
        print("⚠️ No BLAST files found.")
        return None, []

    all_hits = []
    for file in blast_files:
        try:
            df = pd.read_csv(file, sep="\t", header=None, names=[
                "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
            ])
            if df.empty:
                continue
            df["sample"] = os.path.basename(file).split("_blast")[0]
            all_hits.append(df)
        except Exception as e:
            print(f"❌ Error reading {file}: {e}")

    if not all_hits:
        print("⚠️ No valid BLAST data.")
        return None, []

    blast_df = pd.concat(all_hits, ignore_index=True)
    blast_df = blast_df[pd.to_numeric(blast_df["pident"], errors="coerce") >= min_identity]
    blast_df = blast_df[blast_df["stitle"].str.contains("Human herpesvirus 5|Cytomegalovirus Merlin|NC_006273.2", case=False, na=False)]

    rows = []
    for _, row in blast_df.iterrows():
        start = int(min(row["sstart"], row["send"]))
        end = int(max(row["sstart"], row["send"]))
        # include at least one bin
        b0 = start - (start % bin_size)
        bins = range(b0, max(b0 + bin_size, end), bin_size)
        for b in bins:
            rows.append([row["sample"], b, row["pident"]])

    if not rows:
        print("⚠️ No positions to plot.")
        return None, []

    heat_df = pd.DataFrame(rows, columns=["sample", "position", "pident"])
    pivot = heat_df.groupby(["sample", "position"])["pident"].mean().unstack(fill_value=0)

    plt.figure(figsize=(20, 6))
    ax = sns.heatmap(pivot, cmap="coolwarm", cbar_kws={"label": "% Identity"}, vmin=min_identity, vmax=100)
    plt.title("BLAST % Identity Across Genome (100bp bins)")
    plt.xlabel("Genomic Position (binned)")
    plt.ylabel("Sample")

    for region in CONSERVE_REGIONS:
        x_start = region["start"] - (region["start"] % bin_size)
        x_end = region["end"] - (region["end"] % bin_size)
        width = max(bin_size, x_end - x_start)
        rect = patches.Rectangle((x_start, 0), width, len(pivot.index), linewidth=1.5,
                                 edgecolor='black', facecolor='none', clip_on=False)
        ax.add_patch(rect)
        ax.text((x_start + x_end) / 2, -0.5, region["name"], ha="center", va="bottom",
                fontsize=8, color="black", rotation=90)

    os.makedirs(os.path.join(out_dir, "blast_regions"), exist_ok=True)
    full_path = os.path.join(out_dir, "blast_identity_heatmap.png")
    plt.tight_layout()
    plt.savefig(full_path)
    plt.close()

    region_paths = []
    for region in CONSERVE_REGIONS:
        region_df = heat_df[
            (heat_df["position"] >= region["start"]) &
            (heat_df["position"] <= region["end"])
        ]
        if region_df.empty:
            continue
        pivot = region_df.groupby(["sample", "position"])["pident"].mean().unstack(fill_value=0)

        plt.figure(figsize=(12, 4))
        sns.heatmap(pivot, cmap="coolwarm", cbar_kws={"label": "% Identity"}, vmin=min_identity, vmax=100)
        plt.title(f"BLAST % Identity: {region['name']} ({region['start']}-{region['end']})")
        plt.xlabel("Genomic Position (binned)")
        plt.ylabel("Sample")
        out_path = os.path.join(out_dir, "blast_regions", f"{region['name']}_blast_heatmap.png")
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close()
        region_paths.append((region["name"], f"blast_regions/{region['name']}_blast_heatmap.png"))

    return full_path, region_paths

# --- Per-base coverage heatmaps ---
def plot_per_base_depth_heatmaps(per_base_df, out_dir):
    heatmap_dir = os.path.join(out_dir, "per_base_heatmaps")
    os.makedirs(heatmap_dir, exist_ok=True)
    region_map = {r["name"]: (r["start"], r["end"]) for r in CONSERVE_REGIONS}

    for gene, (start, end) in region_map.items():
        data = per_base_df[(per_base_df["gene"] == gene) &
                           (per_base_df["position"] >= start) &
                           (per_base_df["position"] <= end)].copy()
        if data.empty:
            continue
        data["relative_position"] = data["position"] - start + 1
        pivot = data.pivot(index="sample", columns="relative_position", values="depth")

        plt.figure(figsize=(18, 6))
        sns.heatmap(pivot, cmap="viridis", cbar_kws={"label": "Depth"}, vmin=0, vmax=200)
        plt.title(f"Per-base Depth: {gene}")
        plt.xlabel("Position (relative to gene start)")
        plt.ylabel("Sample")

        out_path = os.path.join(heatmap_dir, f"{gene}_per_base_depth.png")
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close()

# --- Genome-wide coverage plots ---
def plot_depth_across_genome(per_base_df, out_dir):
    coverage_dir = os.path.join(out_dir, "genome_coverage")
    os.makedirs(coverage_dir, exist_ok=True)

    for sample in per_base_df["sample"].dropna().unique():
        df = per_base_df[per_base_df["sample"] == sample].sort_values("position")
        if df.empty:
            continue

        plt.figure(figsize=(20, 5))
        plt.plot(df["position"], df["depth"], lw=0.7, color="steelblue")
        plt.yscale("log")

        for region in CONSERVE_REGIONS:
            plt.axvspan(region["start"], region["end"], color="orange", alpha=0.3)
            plt.text((region["start"] + region["end"]) / 2, max(1, df["depth"].max()) * 1.1,
                     region["name"], ha="center", va="bottom", fontsize=9, color="darkred")

        plt.title(f"Per-base Depth Across Genome: {sample}")
        plt.xlabel("Genomic Position")
        plt.ylabel("Depth (log scale)")
        plt.xlim(0, 235000)
        plt.tight_layout()
        plt.savefig(os.path.join(coverage_dir, f"{sample}_genome_coverage.png"))
        plt.close()

# --- Final report assembly ---
def generate_final_cmv_report(gene_df, blast_df, bam_df, per_base_df, fastqc_dir, blast_dir, out_dir):
    """
    API unchanged. Internally, we auto-load and merge all sample CSVs from out_dir
    so figures/tables include *all* samples, not just the last one passed in.
    """
    os.makedirs(out_dir, exist_ok=True)

    # NEW: ensure we use merged data across all samples on disk
    gene_df, blast_df, bam_df, per_base_df = _load_all_metrics_from_disk(
        out_dir=out_dir,
        blast_dir=blast_dir,
        passed_gene_df=gene_df,
        passed_blast_df=blast_df,
        passed_bam_df=bam_df,
        passed_per_base_df=per_base_df
    )

    # --- Gene-level coverage plot ---
    avg_depth_path = None
    if not gene_df.empty:
        plt.figure(figsize=(14, 6))
        sns.barplot(data=gene_df, x="gene", y="avg_depth", hue="sample")
        plt.title("Average Depth per CMV Gene")
        plt.ylabel("Average Depth of Coverage")
        avg_depth_path = os.path.join(out_dir, "01_avg_depth_barplot.png")
        plt.tight_layout()
        plt.savefig(avg_depth_path)
        plt.close()
    else:
        print("⚠️ gene_df is empty: skipping Average Depth plot.")

    # --- Top BLAST hits barplot ---
    blast_bar_path = None
    if not blast_df.empty:
        blast_df["total_hits"] = pd.to_numeric(blast_df["total_hits"], errors="coerce")
        plt.figure(figsize=(14, 6))
        sns.barplot(data=blast_df, x="sample", y="total_hits", hue="species")
        plt.title("Top BLAST Species per Sample and Read")
        plt.xticks(rotation=45, ha="right")
        blast_bar_path = os.path.join(out_dir, "02_top_blast_hits.png")
        plt.tight_layout()
        plt.savefig(blast_bar_path)
        plt.close()
    else:
        print("⚠️ blast_df is empty: skipping Top BLAST hits plot.")

    # --- Other visualizations (BLAST heatmaps, per-base heatmaps, genome coverage) ---
    blast_heatmap_path, region_heatmaps = plot_blast_identity_heatmaps(blast_dir, out_dir)
    if not per_base_df.empty:
        plot_per_base_depth_heatmaps(per_base_df, out_dir)
        plot_depth_across_genome(per_base_df, out_dir)
    else:
        print("⚠️ per_base_df is empty: skipping per-base and genome-wide coverage plots.")

    # Load images for HTML embedding
    per_base_dir = os.path.join(out_dir, "per_base_heatmaps")
    per_base_imgs = sorted(os.listdir(per_base_dir)) if os.path.isdir(per_base_dir) else []
    per_base_tags = "\n".join(f'<h3>{img}</h3><img src="per_base_heatmaps/{img}">' for img in per_base_imgs)

    genome_dir = os.path.join(out_dir, "genome_coverage")
    genome_imgs = sorted(os.listdir(genome_dir)) if os.path.isdir(genome_dir) else []
    genome_tags = "\n".join(f'<h3>{img}</h3><img src="genome_coverage/{img}">' for img in genome_imgs)

    fastqc_htmls = glob.glob(os.path.join(fastqc_dir, "*_fastqc.html"))
    fastqc_links = "\n".join(
        f'<li><a href="{os.path.relpath(f, out_dir)}" target="_blank">{os.path.basename(f)}</a></li>'
        for f in sorted(fastqc_htmls)
    )

    region_heatmap_tags = "\n".join(
        f'<h4>{name}</h4><img src="{path}">' for name, path in region_heatmaps
    )

    # Conditional tags
    avg_depth_img_tag = f'<img src="{os.path.basename(avg_depth_path)}">' if avg_depth_path else "<p>No Average Depth plot generated.</p>"
    blast_bar_img_tag = f'<img src="{os.path.basename(blast_bar_path)}">' if blast_bar_path else "<p>No Top BLAST hits plot generated.</p>"
    full_blast_heatmap_tag = f"<h3>BLAST % Identity Heatmap (Full)</h3><img src='{os.path.basename(blast_heatmap_path)}'>" if blast_heatmap_path else "<p>No full BLAST heatmap.</p>"
    region_block = region_heatmap_tags if region_heatmaps else "<p>No region-specific heatmaps generated.</p>"

    html_path = os.path.join(out_dir, "final_report.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>CMV Diagnostic Report</title>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery.tablesorter/2.31.3/js/jquery.tablesorter.min.js"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/jquery.tablesorter/2.31.3/css/theme.default.min.css">
    <script>$(function(){{ $("table").tablesorter(); }});</script>
    <style>
        body {{ font-family: Arial; padding: 20px; }}
        img {{ max-width: 100%; margin-bottom: 20px; border: 1px solid #ccc; }}
        table {{ margin-bottom: 30px; }}
        th {{ cursor: pointer; }}
        details {{ margin-bottom: 20px; }}
        summary h2 {{ display: inline; }}
    </style>
</head>
<body>
    <h1>📋 CMV Diagnostic Report</h1>

    <details open><summary><h2>1. Average Depth per Gene</h2></summary>
        {avg_depth_img_tag}
        <div class="sortable">{gene_df.to_html(index=False, classes='tablesorter') if not gene_df.empty else "<p>No gene-level data.</p>"}</div>
    </details>

    <details open><summary><h2>2. BLAST Summary</h2></summary>
        {blast_bar_img_tag}
        <div class="sortable">{blast_df.to_html(index=False, classes='tablesorter') if not blast_df.empty else "<p>No BLAST summary data.</p>"}</div>
        {full_blast_heatmap_tag}
        <details><summary><h3>BLAST % Identity by Conserved Region</h3></summary>
            {region_block}
        </details>
    </details>

    <details><summary><h2>3. Per-base Coverage Heatmaps</h2></summary>
        {per_base_tags if per_base_tags else "<p>No per-base heatmaps generated.</p>"}
    </details>

    <details><summary><h2>4. Genome-wide Coverage</h2></summary>
        {genome_tags if genome_tags else "<p>No genome-wide coverage plots generated.</p>"}
    </details>

    <details><summary><h2>5. FastQC Reports</h2></summary>
        <ul>{fastqc_links if fastqc_links else "<p>No FastQC reports found.</p>"}</ul>
    </details>
</body>
</html>""")
    print(f"✅ Final report generated: {html_path}")
