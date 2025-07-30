import os
import re
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

# ---------------------- Loaders (aggregate ALL samples) ----------------------

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
        for col in ("avg_depth", "breadth", "depth", "position", "total_reads"):
            if col in out.columns:
                out[col] = pd.to_numeric(out[col], errors="coerce")
        return out
    return pd.DataFrame()

def _load_all_from_dirs(csv_dir, blast_dir):
    gene_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_gene.csv")))
    per_base_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_per_base.csv")))
    bam_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_summary.csv")))

    # BLAST top-5 (already per sample from summarize_blast_hits)
    blast_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_blast_top5.csv")))
    # blast heatmaps read *_blast.tsv directly in plot function
    def _uniq(col, df):
        return list(df[col].dropna().unique()) if col in df.columns and not df.empty else []
    print(f"ℹ️ gene_df samples: {_uniq('sample', gene_df)}")
    print(f"ℹ️ per_base_df samples: {_uniq('sample', per_base_df)}")
    print(f"ℹ️ bam_df samples: {_uniq('sample', bam_df)}")
    print(f"ℹ️ blast_df samples: {_uniq('sample', blast_df)}")
    return gene_df, blast_df, bam_df, per_base_df

# ---------------------- Alignment Stats (explicit directory) ----------------------

def _first_int(s):
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else 0

def _first_float(s):
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    return float(m.group(1)) if m else 0.0

def _parse_flagstat(path):
    sample = os.path.basename(path).replace("_flagstat.txt", "")
    total = mapped = properly = singletons = duplicates = paired = 0
    mapped_pct = properly_pct = singletons_pct = 0.0
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            L = line.strip()
            if "in total" in L:
                total = _first_int(L)
            elif "duplicates" in L:
                duplicates = _first_int(L)
            elif "mapped (" in L:
                mapped = _first_int(L)
                m = re.search(r"\(([\d\.]+)%", L)
                if m:
                    mapped_pct = float(m.group(1))
            elif "paired in sequencing" in L:
                paired = _first_int(L)
            elif "properly paired" in L:
                properly = _first_int(L)
                m = re.search(r"\(([\d\.]+)%", L)
                if m:
                    properly_pct = float(m.group(1))
            elif "singletons" in L:
                singletons = _first_int(L)
                m = re.search(r"\(([\d\.]+)%", L)
                if m:
                    singletons_pct = float(m.group(1))
    unmapped = max(0, total - mapped)
    dup_pct = (duplicates / total * 100.0) if total else 0.0
    return {
        "sample": sample,
        "total_reads": total,
        "mapped_reads": mapped,
        "unmapped_reads": unmapped,
        "mapped_pct": round(mapped_pct, 2),
        "properly_paired_reads": properly,
        "properly_paired_pct": round(properly_pct, 2),
        "paired_in_sequencing": paired,
        "singletons": singletons,
        "singletons_pct": round(singletons_pct, 2),
        "duplicates": duplicates,
        "duplicates_pct": round(dup_pct, 2),
    }

def _parse_stats_sn(path):
    sample = os.path.basename(path).replace("_stats.txt", "")
    ins_mean = None
    ins_sd = None
    read_len = None
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("SN"):
                continue
            if "insert size average:" in line:
                ins_mean = _first_float(line)
            elif "insert size standard deviation:" in line:
                ins_sd = _first_float(line)
            elif "average length:" in line:
                read_len = _first_float(line)
    return {"sample": sample, "insert_size_mean": ins_mean, "insert_size_sd": ins_sd, "read_length_mean": read_len}

def _load_alignment_stats_from_dir(bam_stats_dir):
    if not bam_stats_dir or not os.path.isdir(bam_stats_dir):
        print("⚠️ BAM Stats directory not provided or does not exist.")
        return pd.DataFrame(), pd.DataFrame()
    flag_rows, sn_rows = [], []
    for f in glob.glob(os.path.join(bam_stats_dir, "*_flagstat.txt")):
        try:
            flag_rows.append(_parse_flagstat(f))
        except Exception as e:
            print(f"❌ Failed to parse flagstat {f}: {e}")
    for s in glob.glob(os.path.join(bam_stats_dir, "*_stats.txt")):
        try:
            sn_rows.append(_parse_stats_sn(s))
        except Exception as e:
            print(f"❌ Failed to parse stats {s}: {e}")
    flag_df = pd.DataFrame(flag_rows) if flag_rows else pd.DataFrame()
    sn_df = pd.DataFrame(sn_rows) if sn_rows else pd.DataFrame()
    if not flag_df.empty and not sn_df.empty:
        flag_df = flag_df.merge(sn_df, on="sample", how="left")
    if "sample" in flag_df.columns:
        flag_df = flag_df.sort_values("sample").reset_index(drop=True)
    return flag_df, sn_df

# ---------------------- Plotting helpers ----------------------

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

def _plot_alignment_figures(align_df, out_dir):
    img = {}
    if align_df.empty:
        return img
    fig_dir = os.path.join(out_dir, "alignment_stats")
    os.makedirs(fig_dir, exist_ok=True)

    # 1) Stacked mapped vs unmapped
    try:
        dfc = align_df[["sample", "mapped_reads", "unmapped_reads"]].set_index("sample")
        order = dfc.index.tolist()
        mapped = dfc["mapped_reads"].values
        unmapped = dfc["unmapped_reads"].values
        plt.figure(figsize=(12, max(3, 0.5 * len(order) + 1)))
        left = [0]*len(order)
        plt.barh(order, mapped, label="Mapped reads")
        plt.barh(order, unmapped, left=mapped, label="Unmapped reads")
        plt.xlabel("Reads")
        plt.title("Mapped vs Unmapped Reads per Sample")
        plt.legend()
        path = os.path.join(fig_dir, "mapped_unmapped_stacked.png")
        plt.tight_layout()
        plt.savefig(path)
        plt.close()
        img["stacked"] = os.path.relpath(path, out_dir)
    except Exception as e:
        print(f"⚠️ Stacked plot failed: {e}")

    # 2) Mapped %
    try:
        plt.figure(figsize=(10, 5))
        sns.barplot(data=align_df, x="sample", y="mapped_pct")
        plt.ylabel("Mapped (%)")
        plt.title("Alignment Rate (Mapped %)")
        plt.xticks(rotation=45, ha="right")
        path = os.path.join(fig_dir, "mapped_percent.png")
        plt.tight_layout()
        plt.savefig(path)
        plt.close()
        img["mapped_pct"] = os.path.relpath(path, out_dir)
    except Exception as e:
        print(f"⚠️ Mapped% plot failed: {e}")

    # 3) Properly paired %
    try:
        plt.figure(figsize=(10, 5))
        sns.barplot(data=align_df, x="sample", y="properly_paired_pct")
        plt.ylabel("Properly paired (%)")
        plt.title("Properly Paired Reads (%)")
        plt.xticks(rotation=45, ha="right")
        path = os.path.join(fig_dir, "properly_paired_percent.png")
        plt.tight_layout()
        plt.savefig(path)
        plt.close()
        img["properly_paired_pct"] = os.path.relpath(path, out_dir)
    except Exception as e:
        print(f"⚠️ Properly paired% plot failed: {e}")

    # 4) Insert size mean ± SD (if available)
    if "insert_size_mean" in align_df.columns and align_df["insert_size_mean"].notna().any():
        try:
            df_is = align_df[["sample", "insert_size_mean", "insert_size_sd"]].copy()
            plt.figure(figsize=(10, 5))
            ax = sns.barplot(data=df_is, x="sample", y="insert_size_mean")
            yerr = df_is["insert_size_sd"].fillna(0.0).values
            for i, patch in enumerate(ax.patches):
                x = patch.get_x() + patch.get_width()/2
                y = patch.get_height()
                ax.errorbar(x, y, yerr=yerr[i] if i < len(yerr) else 0.0, fmt='none', capsize=3)
            plt.ylabel("Insert size (mean ± SD)")
            plt.title("Insert Size by Sample")
            plt.xticks(rotation=45, ha="right")
            path = os.path.join(fig_dir, "insert_size.png")
            plt.tight_layout()
            plt.savefig(path)
            plt.close()
            img["insert_size"] = os.path.relpath(path, out_dir)
        except Exception as e:
            print(f"⚠️ Insert size plot failed: {e}")
    return img

# ---------------------- Final Report (NEW SIGNATURE) ----------------------

def generate_final_cmv_report(csv_dir, fastqc_dir, blast_dir, out_dir, bam_stats_dir):
    """
    Build the CMV report using directories so ALL samples are included.
    Section numbering:
      1. Alignment Statistics (BWA-MEM)
      2. Average Depth per Gene
      3. BLAST Summary
      4. Per-base Coverage Heatmaps
      5. Genome-wide Coverage
      6. FastQC Reports
    """
    os.makedirs(out_dir, exist_ok=True)

    # Load all CSVs across samples
    gene_df, blast_df, bam_df, per_base_df = _load_all_from_dirs(csv_dir, blast_dir)

    # Load alignment stats explicitly from user-provided BAM_Stats directory
    align_df, _ = _load_alignment_stats_from_dir(bam_stats_dir)
    align_imgs = _plot_alignment_figures(align_df, out_dir)

    # 2) Gene-level coverage plot
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

    # 3) Top BLAST hits plot
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

    # BLAST identity heatmaps
    blast_heatmap_path, region_heatmaps = plot_blast_identity_heatmaps(blast_dir, out_dir)

    # 4 & 5) Coverage visuals
    if not per_base_df.empty:
        plot_per_base_depth_heatmaps(per_base_df, out_dir)
        plot_depth_across_genome(per_base_df, out_dir)
    else:
        print("⚠️ per_base_df is empty: skipping per-base and genome-wide coverage plots.")

    # Collect image tags
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

    # Alignment block (section 1)
    align_table_html = (align_df.rename(columns={
        "mapped_pct":"mapped_%", "properly_paired_pct":"properly_paired_%", "duplicates_pct":"duplicates_%", "singletons_pct":"singletons_%"
    }).to_html(index=False, classes='tablesorter')) if not align_df.empty else "<p>No alignment stats found.</p>"

    align_imgs_tags = []
    if "stacked" in align_imgs:
        align_imgs_tags.append(f'<h3>Mapped vs Unmapped</h3><img src="{align_imgs["stacked"]}">')
    if "mapped_pct" in align_imgs:
        align_imgs_tags.append(f'<h3>Alignment Rate (Mapped %)</h3><img src="{align_imgs["mapped_pct"]}">')
    if "properly_paired_pct" in align_imgs:
        align_imgs_tags.append(f'<h3>Properly Paired (%)</h3><img src="{align_imgs["properly_paired_pct"]}">')
    if "insert_size" in align_imgs:
        align_imgs_tags.append(f'<h3>Insert Size (mean ± SD)</h3><img src="{align_imgs["insert_size"]}">')
    align_imgs_block = "\n".join(align_imgs_tags) if align_imgs_tags else "<p>No alignment figures generated.</p>"

    # HTML assembly with requested numbering
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

    <details open><summary><h2>1. Alignment Statistics (BWA‑MEM)</h2></summary>
        <div class="sortable">{align_table_html}</div>
        {align_imgs_block}
    </details>

    <details open><summary><h2>2. Average Depth per Gene</h2></summary>
        {avg_depth_img_tag}
        <div class="sortable">{gene_df.to_html(index=False, classes='tablesorter') if not gene_df.empty else "<p>No gene-level data.</p>"}</div>
    </details>

    <details open><summary><h2>3. BLAST Summary</h2></summary>
        {blast_bar_img_tag}
        <div class="sortable">{blast_df.to_html(index=False, classes='tablesorter') if not blast_df.empty else "<p>No BLAST summary data.</p>"}</div>
        {full_blast_heatmap_tag}
        <details><summary><h3>BLAST % Identity by Conserved Region</h3></summary>
            {region_block}
        </details>
    </details>

    <details><summary><h2>4. Per-base Coverage Heatmaps</h2></summary>
        {per_base_tags if per_base_tags else "<p>No per-base heatmaps generated.</p>"}
    </details>

    <details><summary><h2>5. Genome-wide Coverage</h2></summary>
        {genome_tags if genome_tags else "<p>No genome-wide coverage plots generated.</p>"}
    </details>

    <details><summary><h2>6. FastQC Reports</h2></summary>
        <ul>{fastqc_links if fastqc_links else "<p>No FastQC reports found.</p>"}</ul>
    </details>
</body>
</html>""")
    print(f"✅ Final report generated: {html_path}")
