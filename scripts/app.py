# app.py
# Streamlit app to interactively explore CMV results for ONE OR MORE canonical samples ("bases")

import os
import re
import io
import glob
import gzip
import zipfile
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import streamlit.components.v1 as components

# --- optional: clear Streamlit caches on startup if env flag is set ---
# Clear Streamlit caches on startup when launched from the GUI
if os.environ.get("CLEAR_CACHE_ON_START") == "1":
    try:
        st.cache_data.clear()
        st.cache_resource.clear()
    except Exception:
        pass


# ----------------------- Conserved regions -----------------------
CONSERVE_REGIONS = [
    {"name": "IE1",   "start": 172328, "end": 174090},
    {"name": "UL54",  "start": 78193,  "end": 81922},
    {"name": "gpUL55","start": 82065,  "end": 84789},
    {"name": "UL75",  "start": 109223, "end": 111452},
    {"name": "UL97",  "start": 141797, "end": 143921},
    {"name": "UL83",  "start": 120655, "end": 122341},
]

# ----------------------- Sample canonicalization -----------------------
_SUFFIX_PATTERNS = [
    r'(?:_R?[12](?:_00[1-9])?)$',          # _R1, _R2, _R1_001, ...
    r'(?:_[12])$',                          # _1, _2
    r'(?:_sorted(?:_sorted)*)$',            # _sorted, _sorted_sorted, ...
    r'(?:_sort(?:ed)?)$',                   # _sort/_sorted
    r'(?:_dedup(?:ed)?)$',                  # _dedup/_deduped
    r'(?:_mark(?:ed)?dups?)$',              # _markdup/_markeddups
    r'(?:_rmdup[s]?)$',                     # _rmdup
    r'(?:_trim(?:med)?(?:_paired)?)$',      # _trim/_trimmed/_trimmed_paired
    r'(?:_filtered)$',                      # _filtered
    r'(?:\.bam)$',                          # .bam
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


def add_base_column(df: pd.DataFrame) -> pd.DataFrame:
    if df is not None and not df.empty and "sample" in df.columns:
        df = df.copy()
        df["base"] = df["sample"].astype(str).apply(canonical_sample)
    return df

# ----------------------- Data loaders -----------------------
def _read_csv_list(files, **read_kwargs) -> pd.DataFrame:
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, **read_kwargs)
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            st.warning(f"Error reading {f}: {e}")
    if not dfs:
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)
    for col in ("avg_depth", "breadth", "depth", "position", "total_reads", "total_hits"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


@st.cache_data(show_spinner=False)
def load_all(csv_dir: str, blast_dir: str):
    gene_df     = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_gene.csv")))
    per_base_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_per_base.csv")))
    bam_df      = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_summary.csv")))
    blast_top5  = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_blast_top5.csv")))

    # BLAST hits (per-sample *_blast.tsv)
    blast_files = glob.glob(os.path.join(blast_dir, "*_blast.tsv"))
    hits = []
    for file in blast_files:
        try:
            df = pd.read_csv(file, sep="\t", header=None, names=[
                "qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore","stitle"
            ])
            if df.empty:
                continue
            df["sample"] = os.path.basename(file).split("_blast")[0]
            df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
            hits.append(df)
        except Exception as e:
            st.warning(f"Error reading {file}: {e}")
    blast_hits = pd.concat(hits, ignore_index=True) if hits else pd.DataFrame()

    # attach canonical base IDs
    gene_df     = add_base_column(gene_df)
    per_base_df = add_base_column(per_base_df)
    bam_df      = add_base_column(bam_df)
    blast_top5  = add_base_column(blast_top5)
    blast_hits  = add_base_column(blast_hits)
    return gene_df, per_base_df, bam_df, blast_top5, blast_hits


@st.cache_data(show_spinner=False)
def load_alignment_stats(bam_stats_dir: str) -> pd.DataFrame:
    def _first_int(s):
        m = re.search(r"(\d+)", s);  return int(m.group(1)) if m else 0
    def _first_float(s):
        m = re.search(r"(\d+(?:\.\d+)?)", s);  return float(m.group(1)) if m else 0.0

    flag_rows, sn_rows = [], []
    for f in glob.glob(os.path.join(bam_stats_dir, "*_flagstat.txt")):
        try:
            sample = os.path.basename(f).replace("_flagstat.txt", "")
            total = mapped = properly = singletons = duplicates = paired = 0
            mapped_pct = properly_pct = singletons_pct = 0.0
            with open(f, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    L = line.strip()
                    if "in total" in L:
                        total = _first_int(L)
                    elif "duplicates" in L:
                        duplicates = _first_int(L)
                    elif "mapped (" in L:
                        mapped = _first_int(L)
                        m = re.search(r"\(([\d\.]+)%", L)
                        if m: mapped_pct = float(m.group(1))
                    elif "paired in sequencing" in L:
                        paired = _first_int(L)
                    elif "properly paired" in L:
                        properly = _first_int(L)
                        m = re.search(r"\(([\d\.]+)%", L)
                        if m: properly_pct = float(m.group(1))
                    elif "singletons" in L:
                        singletons = _first_int(L)
                        m = re.search(r"\(([\d\.]+)%", L)
                        if m: singletons_pct = float(m.group(1))
            unmapped = max(0, total - mapped)
            dup_pct = (duplicates / total * 100.0) if total else 0.0
            flag_rows.append({
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
            })
        except Exception as e:
            st.warning(f"Failed to parse flagstat {f}: {e}")

    for s in glob.glob(os.path.join(bam_stats_dir, "*_stats.txt")):
        try:
            sample = os.path.basename(s).replace("_stats.txt", "")
            ins_mean = ins_sd = read_len = None
            with open(s, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith("SN"): continue
                    if "insert size average:" in line:              ins_mean = _first_float(line)
                    elif "insert size standard deviation:" in line: ins_sd   = _first_float(line)
                    elif "average length:" in line:                 read_len = _first_float(line)
            sn_rows.append({"sample": sample, "insert_size_mean": ins_mean, "insert_size_sd": ins_sd, "read_length_mean": read_len})
        except Exception as e:
            st.warning(f"Failed to parse stats {s}: {e}")

    flag_df = pd.DataFrame(flag_rows)
    sn_df   = pd.DataFrame(sn_rows)
    if not flag_df.empty and not sn_df.empty:
        flag_df = flag_df.merge(sn_df, on="sample", how="left")
    if "sample" in flag_df.columns:
        flag_df = flag_df.sort_values("sample").reset_index(drop=True)

    # attach canonical base IDs
    flag_df = add_base_column(flag_df)
    return flag_df


# ----------------------- NEW: VAF loaders & helpers -----------------------
@st.cache_data(show_spinner=False)
def load_vaf_tables(vaf_dir: str) -> pd.DataFrame:
    """Load LoFreq VAF tables: expects *_vaf.tsv with header CHROM POS REF ALT DP AD_ALT AF."""
    if not vaf_dir or not os.path.isdir(vaf_dir):
        return pd.DataFrame()
    paths = glob.glob(os.path.join(vaf_dir, "*_vaf.tsv"))
    rows = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t")
            if df.empty:
                continue
            sample = os.path.basename(p).replace("_vaf.tsv", "")
            df["sample"] = sample
            rows.append(df)
        except Exception as e:
            st.warning(f"Error reading VAF file {p}: {e}")
    if not rows:
        return pd.DataFrame()
    vaf = pd.concat(rows, ignore_index=True)
    # Numeric coercion
    for c in ("POS", "DP", "AD_ALT", "AF"):
        if c in vaf.columns:
            vaf[c] = pd.to_numeric(vaf[c], errors="coerce")
    vaf = add_base_column(vaf)
    return vaf


def infer_mixed_strains(vaf_df: pd.DataFrame) -> str:
    """Very simple heuristic: mixed strains if there are enough mid-VAF variants.
    Returns one of: 'Likely', 'Possible', 'Unlikely'."""
    if vaf_df is None or vaf_df.empty:
        return "Unknown"
    df = vaf_df.dropna(subset=["AF"]).copy()
    n_total = len(df)
    n_mid = int(((df["AF"] >= 0.2) & (df["AF"] <= 0.8)).sum())
    if n_total == 0:
        return "Unknown"
    frac = n_mid / n_total
    if n_mid >= 10 and frac >= 0.15:
        return "Likely"
    if n_mid >= 3 and frac >= 0.05:
        return "Possible"
    return "Unlikely"


def detect_non_cmv_species(blast_hits_all: pd.DataFrame, base: str) -> bool:
    if blast_hits_all is None or blast_hits_all.empty or "base" not in blast_hits_all.columns:
        return False
    sub = blast_hits_all[blast_hits_all["base"] == base]
    if sub.empty:
        return False
    # Look for non-CMV strings in stitle
    cmv_pat = re.compile(r"Human\s+herpesvirus\s*5|\bHCMV\b|Cytomegalovirus|NC_006273\.2", re.IGNORECASE)
    non = sub[~sub["stitle"].fillna("").str.contains(cmv_pat)]
    # If there are a reasonable number of non-CMV hits, flag
    return len(non) >= 50  # threshold can be adjusted


def plot_vaf_by_position(vaf_df: pd.DataFrame):
    if vaf_df is None or vaf_df.empty:
        return None
    fig = go.Figure()
    for sample, sdf in vaf_df.groupby("sample"):
        sdf = sdf.sort_values("POS")
        fig.add_trace(go.Scatter(x=sdf["POS"], y=sdf["AF"], mode="markers",
                                 name=sample, hovertemplate="POS=%{x}<br>AF=%{y:.3f}<extra>"))
    fig.update_layout(title="Variant Allele Frequency (AF) across genome",
                      xaxis_title="Position (bp)", yaxis_title="AF (0-1)",
                      yaxis=dict(range=[0,1]))
    return fig


def plot_af_histogram(vaf_df: pd.DataFrame):
    if vaf_df is None or vaf_df.empty:
        return None
    fig = px.histogram(vaf_df, x="AF", color="sample", nbins=40, title="Distribution of AF values")
    fig.update_xaxes(range=[0,1])
    return fig


# ----------------------- Filters & helpers -----------------------
def filter_blast_merlin(blast_hits: pd.DataFrame, min_identity=85.0) -> pd.DataFrame:
    if blast_hits.empty: return blast_hits
    df = blast_hits.copy()
    df = df[pd.to_numeric(df["pident"], errors="coerce").notna()]
    title_mask = df["stitle"].str.contains(
        "Human herpesvirus 5|Human cytomegalovirus|Cytomegalovirus Merlin|NC_006273.2",
        case=False, na=False
    )
    df = df[title_mask]
    return df[df["pident"] >= float(min_identity)]


def filter_by_bases(df: pd.DataFrame, bases_selected: list[str]) -> pd.DataFrame:
    if df is None or df.empty or "base" not in df.columns or not bases_selected:
        return df
    return df[df["base"].isin(bases_selected)].copy()


# ----------------------- Plotting -----------------------
def plot_alignment_stacked(df: pd.DataFrame):
    if df.empty: return None
    fig = go.Figure()
    fig.add_bar(y=df["sample"], x=df["mapped_reads"], orientation="h", name="Mapped reads")
    fig.add_bar(y=df["sample"], x=df["unmapped_reads"], orientation="h", name="Unmapped reads")
    fig.update_layout(barmode="stack", xaxis_title="Reads", yaxis_title="Sample", title="Mapped vs Unmapped Reads per Sample")
    return fig


def plot_alignment_pct(df: pd.DataFrame, col: str, title: str):
    if df.empty or col not in df: return None
    fig = px.bar(df, x="sample", y=col, title=title)
    fig.update_yaxes(title="%"); fig.update_xaxes(tickangle=45)
    return fig


def plot_avg_depth(gene_df: pd.DataFrame):
    if gene_df.empty: return None
    fig = px.bar(gene_df, x="gene", y="avg_depth", color="sample", title="Average Depth per CMV Gene")
    fig.update_yaxes(title="Average Depth"); fig.update_xaxes(title="Gene")
    return fig


def plot_blast_full_heatmap(blast_hits: pd.DataFrame, bin_size=100):
    if blast_hits.empty: return None
    rows = []
    for _, r in blast_hits.iterrows():
        s0 = int(min(r["sstart"], r["send"])); s1 = int(max(r["sstart"], r["send"]))
        b0 = s0 - (s0 % bin_size)
        for b in range(b0, max(b0 + bin_size, s1), bin_size):
            rows.append([r["sample"], b, float(r["pident"])])
    if not rows: return None
    heat_df = pd.DataFrame(rows, columns=["sample", "position", "pident"])
    pivot = heat_df.groupby(["sample", "position"]).pident.mean().unstack(fill_value=0).sort_index(axis=1)
    fig = px.imshow(pivot, aspect="auto", color_continuous_scale="RdBu", zmin=0, zmax=100,
                    labels=dict(color="% Identity"), title=f"BLAST % Identity Across Genome ({bin_size} bp bins)")
    # conserved regions as vertical rectangles
    shapes = []
    for reg in CONSERVE_REGIONS:
        x0 = reg["start"] - (reg["start"] % bin_size)
        x1 = reg["end"]   - (reg["end"]   % bin_size)
        shapes.append(dict(type="rect", xref="x", yref="paper", x0=x0, x1=x1, y0=0, y1=1,
                           line=dict(color="black"), fillcolor="rgba(0,0,0,0)"))
    fig.update_layout(shapes=shapes)
    return fig


def plot_blast_region_heatmap(blast_hits: pd.DataFrame, region: dict):
    if blast_hits.empty: return None
    rstart, rend = int(region["start"]), int(region["end"])
    region_len = rend - rstart + 1
    columns = np.arange(rstart, rend + 1)
    samples = sorted(blast_hits["sample"].unique())

    sums = np.zeros((len(samples), region_len), dtype=float)
    counts = np.zeros_like(sums, dtype=int)
    sidx = {s: i for i, s in enumerate(samples)}

    for _, r in blast_hits.iterrows():
        s0 = int(min(r["sstart"], r["send"])); s1 = int(max(r["sstart"], r["send"]))
        left, right = max(s0, rstart), min(s1, rend)
        if left > right: continue
        i = sidx.get(r["sample"])
        a, b = left - rstart, right - rstart + 1
        p = float(r["pident"]) 
        sums[i, a:b] += p
        counts[i, a:b] += 1

    means = np.divide(sums, counts, out=np.zeros_like(sums, dtype=float), where=counts > 0)
    pivot = pd.DataFrame(means, index=samples, columns=columns)
    fig = px.imshow(pivot, aspect="auto", color_continuous_scale="RdBu", zmin=0, zmax=100,
                    labels=dict(color="% Identity"),
                    title=f"BLAST % Identity: {region['name']} ({rstart}-{rend})")
    return fig


def plot_per_base_heatmap(per_base_df: pd.DataFrame, gene: str, bin_size: int | None = None):
    region_map = {r["name"]: (r["start"], r["end"]) for r in CONSERVE_REGIONS}
    start, end = region_map[gene]
    data = per_base_df[(per_base_df["gene"] == gene) & (per_base_df["position"] >= start) & (per_base_df["position"] <= end)].copy()
    if data.empty: return None
    if bin_size and bin_size > 1:
        data["bin_pos"] = (data["position"] // bin_size) * bin_size
        pivot = data.groupby(["sample", "bin_pos"]).depth.mean().unstack(fill_value=0).sort_index(axis=1)
        title = f"Per-base Depth (binned {bin_size} bp): {gene} [{start}-{end}]"
    else:
        pivot = data.pivot(index="sample", columns="position", values="depth").fillna(0).sort_index(axis=1)
        title = f"Per-base Depth: {gene} [{start}-{end}]"
    fig = px.imshow(pivot, aspect="auto", color_continuous_scale="Viridis", labels=dict(color="Depth"), title=title)
    return fig


def plot_genome_coverage(per_base_df: pd.DataFrame):
    if per_base_df.empty: return None
    fig = go.Figure()
    for sample, df in per_base_df.groupby("sample"):
        df = df.sort_values("position")
        fig.add_trace(go.Scatter(x=df["position"], y=df["depth"], mode="lines", name=sample))
    for region in CONSERVE_REGIONS:
        fig.add_vrect(x0=region["start"], x1=region["end"], fillcolor="orange", opacity=0.2, line_width=0,
                      annotation_text=region["name"], annotation_position="top")
    fig.update_layout(title="Per-base Depth Across Genome", xaxis_title="Genomic Position",
                      yaxis_title="Depth (log)", yaxis_type="log", xaxis_range=[0, 235000])
    return fig


# ----------------------- FastQC helpers -----------------------
def parse_fastqc_zip_summary(zip_path: str) -> pd.DataFrame:
    sample = os.path.basename(zip_path).replace("_fastqc.zip", "")
    rows = []
    try:
        with zipfile.ZipFile(zip_path) as zf:
            inner = [n for n in zf.namelist() if n.endswith("/summary.txt")]
            if not inner: return pd.DataFrame()
            with zf.open(inner[0]) as fh:
                for line in io.TextIOWrapper(fh, encoding="utf-8", errors="ignore"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        status, module, _filename = parts[0], parts[1], parts[2]
                        rows.append({"sample": sample, "module": module, "status": status})
    except Exception as e:
        st.warning(f"Failed to read {zip_path}: {e}")
    return pd.DataFrame(rows)


def embed_fastqc_html(path: str):
    try:
        if path.endswith(".zip"):
            with zipfile.ZipFile(path) as zf:
                html_files = [f for f in zf.namelist() if f.endswith("fastqc_report.html")]
                if html_files:
                    with zf.open(html_files[0]) as f:
                        html_content = f.read().decode("utf-8", errors="ignore")
                        components.html(html_content, height=800, scrolling=True)
                        return
        elif path.endswith(".html"):
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                html_content = f.read()
                components.html(html_content, height=800, scrolling=True)
                return
        st.info("Could not find a FastQC HTML to render.")
    except Exception as e:
        st.warning(f"Preview failed for {path}: {e}")


# ----------------------- UI -----------------------
st.set_page_config(page_title="CMV Interactive Report", layout="wide")
st.title("📊 CMV Interactive Report")

with st.sidebar:
    st.header("Inputs")
    csv_dir = st.text_input("CSV directory (gene/per_base/summary/blast_top5)")
    blast_dir = st.text_input("BLAST directory (*_blast.tsv)")
    vaf_dir = st.text_input("VAF directory (*_vaf.tsv)")  # NEW
    bam_stats_dir = st.text_input("BAM stats directory (*_flagstat.txt, *_stats.txt)")
    fastqc_dir = st.text_input("FastQC directory (*_fastqc.html, *_fastqc.zip)")
    min_identity = st.number_input("BLAST min % identity", min_value=0.0, max_value=100.0, value=85.0, step=1.0)
    bin_size = st.number_input("BLAST heatmap bin (bp)", min_value=10, max_value=1000, value=100, step=10)

# autoload when directories are provided
if csv_dir and blast_dir:
    gene_df, per_base_df, bam_df, blast_top5, blast_hits_all = load_all(csv_dir, blast_dir)
    align_df = load_alignment_stats(bam_stats_dir) if bam_stats_dir else pd.DataFrame()
    vaf_df = load_vaf_tables(vaf_dir) if vaf_dir else pd.DataFrame()

    # ---------------- FIX: Build base list from core CMV CSVs only ----------------
    core_sets = []
    for df in (gene_df, per_base_df, bam_df):
        if df is not None and not df.empty and "base" in df.columns:
            core_sets.append(set(df["base"].dropna().unique()))
    core_bases = set.union(*core_sets) if core_sets else set()

    # Fallback to BLAST if no CSVs present (edge case)
    if not core_bases:
        if blast_top5 is not None and not blast_top5.empty and "base" in blast_top5.columns:
            core_bases |= set(blast_top5["base"].dropna().unique())
        if blast_hits_all is not None and not blast_hits_all.empty and "base" in blast_hits_all.columns:
            core_bases |= set(blast_hits_all["base"].dropna().unique())

    all_bases = sorted(core_bases)

    # Filter every dataframe to these bases so plots & tables stay consistent
    def keep_bases(df):
        if df is None or df.empty or "base" not in df.columns:
            return df
        return df[df["base"].isin(core_bases)].copy()

    gene_df        = keep_bases(gene_df)
    per_base_df    = keep_bases(per_base_df)
    bam_df         = keep_bases(bam_df)
    blast_top5     = keep_bases(blast_top5)
    blast_hits_all = keep_bases(blast_hits_all)
    align_df       = keep_bases(align_df)
    vaf_df         = keep_bases(vaf_df)
    # ------------------------------------------------------------------------------

    if not all_bases:
        st.warning("No samples found across inputs. Check your folder paths and file naming.")
        st.stop()

    # (Optional) tiny debug helper
    with st.sidebar.expander("Debug: bases detected by source", expanded=False):
        def names(df):
            return sorted(df["base"].dropna().unique()) if (df is not None and not df.empty and "base" in df.columns) else []
        st.write({
            "gene.csv": names(gene_df),
            "per_base.csv": names(per_base_df),
            "summary.csv": names(bam_df),
            "blast_top5.csv": names(blast_top5),
            "*.tsv (BLAST hits)": names(blast_hits_all),
            "BAM stats": names(align_df),
            "VAF tables": names(vaf_df),
        })

    with st.sidebar:
        bases_selected = st.multiselect(
            "Select sample(s)",
            options=all_bases,
            default=[all_bases[0]] if all_bases else [],
            help="Pick one or more canonical sample IDs. All technical variants for each base are included automatically."
        )

        # show which technical variants will be included
        if bases_selected:
            st.markdown("**Files Included:**")
            for b in bases_selected:
                variants_b = sorted(set(
                    s for df in [gene_df, per_base_df, bam_df, blast_top5, blast_hits_all, align_df, vaf_df]
                    if df is not None and not df.empty
                    for s, bb in zip(df["sample"], df["base"])
                    if bb == b
                ))
                st.caption(f"• {b}\n" + ("\n".join(variants_b) if variants_b else "(no variants detected)"))
        else:
            st.info("Select at least one base sample.")
            st.stop()

    # Filter by chosen base(s)
    f_gene  = filter_by_bases(gene_df, bases_selected)
    f_per   = filter_by_bases(per_base_df, bases_selected)
    f_top5  = filter_by_bases(blast_top5, bases_selected)
    f_align = filter_by_bases(align_df, bases_selected)
    sel_hits = filter_by_bases(blast_hits_all, bases_selected)
    sel_hits = filter_blast_merlin(sel_hits, min_identity=float(min_identity))
    f_vaf   = filter_by_bases(vaf_df, bases_selected)

    # 1) Alignment
    with st.expander("1. Alignment Statistics (BWA-MEM)", expanded=True):
        if not f_align.empty:
            st.dataframe(f_align.drop(columns=["base"]), use_container_width=True)
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                fig = plot_alignment_stacked(f_align)
                if fig: st.plotly_chart(fig, use_container_width=True)
            with c2:
                fig = plot_alignment_pct(f_align, "mapped_pct", "Alignment Rate (Mapped %)")
                if fig: st.plotly_chart(fig, use_container_width=True)
            with c3:
                fig = plot_alignment_pct(f_align, "properly_paired_pct", "Properly Paired Reads (%)")
                if fig: st.plotly_chart(fig, use_container_width=True)
            with c4:
                if "insert_size_mean" in f_align.columns and f_align["insert_size_mean"].notna().any():
                    fig = px.bar(f_align, x="sample", y="insert_size_mean", error_y="insert_size_sd",
                                 title="Insert Size (mean ± SD)")
                    st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No alignment stats found for this sample.")

    # 2) Average Depth per Gene
    with st.expander("2. Average Depth per Gene", expanded=True):
        if not f_gene.empty:
            st.dataframe(f_gene.drop(columns=["base"]), use_container_width=True)
            fig = plot_avg_depth(f_gene)
            if fig: st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No gene-level coverage data for this sample.")

    # 3) BLAST Summary + Heatmaps
    with st.expander("3. BLAST Summary & Identity Heatmaps", expanded=True):
        if not f_top5.empty:
            st.subheader("Top BLAST species per sample")
            st.dataframe(f_top5.drop(columns=["base"]), use_container_width=True)
            if "species" in f_top5.columns and "total_hits" in f_top5.columns:
                fig = px.bar(f_top5, x="sample", y="total_hits", color="species",
                             title="Top BLAST Species per Sample and Read")
                fig.update_xaxes(tickangle=45)
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No BLAST summary (top5) for this sample.")

        if not sel_hits.empty:
            st.subheader("BLAST % identity across genome")
            fig = plot_blast_full_heatmap(sel_hits, bin_size=int(bin_size))
            if fig: st.plotly_chart(fig, use_container_width=True)

            st.subheader("BLAST % identity by conserved region")
            tabs = st.tabs([r["name"] for r in CONSERVE_REGIONS])
            for tab, region in zip(tabs, CONSERVE_REGIONS):
                with tab:
                    fig = plot_blast_region_heatmap(sel_hits, region=region)
                    if fig:
                        st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No BLAST hits for this sample (after CMV/Merlin filter & identity threshold).")

    # 4) Per-base Coverage Heatmaps
    with st.expander("4. Per-base Coverage Heatmaps", expanded=False):
        if not f_per.empty:
            tabs = st.tabs([r["name"] for r in CONSERVE_REGIONS])
            for tab, region in zip(tabs, CONSERVE_REGIONS):
                with tab:
                    fig = plot_per_base_heatmap(f_per, gene=region["name"], bin_size=None)
                    if fig: st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No per-base depth data for this sample.")

    # 5) Genome-wide Coverage
    with st.expander("5. Genome-wide Coverage", expanded=False):
        if not f_per.empty:
            fig = plot_genome_coverage(f_per)
            if fig: st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No per-base depth data for genome-wide view.")

    # 6) FastQC Reports
    with st.expander("6. FastQC Reports", expanded=False):
        if fastqc_dir and os.path.isdir(fastqc_dir):
            selected_bases_set = set(bases_selected)
            def _is_selected_fastqc(path: str) -> bool:
                fname = os.path.basename(path)
                sample = fname.replace("_fastqc.html", "").replace("_fastqc.zip", "")
                return canonical_sample(sample) in selected_bases_set

            html_files_all = sorted(glob.glob(os.path.join(fastqc_dir, "*_fastqc.html")))
            zip_files_all  = sorted(glob.glob(os.path.join(fastqc_dir, "*_fastqc.zip")))
            html_files = [p for p in html_files_all if _is_selected_fastqc(p)]
            zip_files  = [p for p in zip_files_all  if _is_selected_fastqc(p)]

            # Parse module summaries from ZIPs
            all_sum = []
            for z in zip_files:
                dfz = parse_fastqc_zip_summary(z)
                if not dfz.empty:
                    all_sum.append(dfz)
            if all_sum:
                df_sum = pd.concat(all_sum, ignore_index=True)
                df_sum["base"] = df_sum["sample"].apply(canonical_sample)
                df_sum = df_sum[df_sum["base"].isin(bases_selected)].drop(columns=["base"])
                st.subheader("FastQC module summary (from ZIPs)")
                pivot = df_sum.pivot_table(index=["sample"], columns="module", values="status",
                                           aggfunc=lambda x: ", ".join(sorted(set(x))))
                st.dataframe(pivot.fillna("N/A"), use_container_width=True)

            reports = html_files + zip_files
            if reports:
                choice = st.selectbox("Preview a FastQC report", options=["(select)"] + reports)
                if choice and choice != "(select)":
                    st.markdown(f"**Preview:** {os.path.basename(choice)}")
                    embed_fastqc_html(choice)

                st.subheader("Download reports")
                c1, c2 = st.columns(2)
                with c1:
                    for f in html_files:
                        name = os.path.basename(f)
                        with open(f, "rb") as fh:
                            st.download_button(label=f"Download {name}", data=fh.read(),
                                               file_name=name, mime="text/html")
                with c2:
                    for z in zip_files:
                        name = os.path.basename(z)
                        with open(z, "rb") as fh:
                            st.download_button(label=f"Download {name}", data=fh.read(),
                                               file_name=name, mime="application/zip")

                st.caption("Inline preview may miss some styling/assets. Use the download buttons for full-fidelity viewing.")
            else:
                st.info("No FastQC HTML or ZIP files found.")
        else:
            st.info("Provide a valid FastQC directory to list reports.")

    # 7) Variants (LoFreq VAF)
    with st.expander("7. Variants (LoFreq VAF)", expanded=True):
        if f_vaf is None or f_vaf.empty:
            st.info("No VAF tables loaded. Provide a VAF directory containing *_vaf.tsv files.")
        else:
            # Show table
            st.subheader("VAF table (selected samples)")
            st.dataframe(
                f_vaf[[c for c in ["sample","CHROM","POS","REF","ALT","DP","AD_ALT","AF","base"] if c in f_vaf.columns]]
                .sort_values(["sample","POS"]).reset_index(drop=True)
                .drop(columns=["base"]),
                use_container_width=True,
            )

            # Plots
            c1, c2 = st.columns(2)
            with c1:
                fig = plot_vaf_by_position(f_vaf)
                if fig: st.plotly_chart(fig, use_container_width=True)
            with c2:
                fig = plot_af_histogram(f_vaf)
                if fig: st.plotly_chart(fig, use_container_width=True)

            # Co-infection / mixed strains flags per base
            st.subheader("Co-infection / Mixed-strain heuristic")
            for b in bases_selected:
                vf_b = f_vaf[f_vaf["base"] == b] if "base" in f_vaf.columns else pd.DataFrame()
                verdict = infer_mixed_strains(vf_b)
                non_cmv = detect_non_cmv_species(blast_hits_all, b)
                msg = f"**{b}:** Mixed HCMV strains: **{verdict}**"
                if non_cmv:
                    msg += "  ·  Non-CMV viral reads detected in BLAST: **Yes**"
                else:
                    msg += "  ·  Non-CMV viral reads detected in BLAST: **No/Low**"
                if verdict == "Likely" or non_cmv:
                    st.warning(msg)
                elif verdict == "Possible":
                    st.info(msg)
                else:
                    st.success(msg)
else:
    st.info("Enter at least CSV and BLAST directories to load data.")
