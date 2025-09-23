# app.py
# Streamlit app to explore CMV results for one or more canonical samples ("bases")
# - No "breadth" metrics
# - Uses selected samples across ALL sections
# - E-value flags overall and per conserved region (Merlin), with significance comments
# - Summary table: removes "Overall", adds "Mixed strains"
# - No error banner for "No CMV hits"

import os
import re
import io
import glob
import zipfile
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import streamlit.components.v1 as components

# ----------------------- Folder picker & state sync -----------------------
def _pick_directory_dialog(initial: str | None = None) -> str | None:
    """Open a native folder picker (tkinter). Returns path or None on cancel/failure."""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        try:
            root.call('wm', 'attributes', '.', '-topmost', True)
        except Exception:
            pass
        directory = filedialog.askdirectory(initialdir=initial or os.getcwd())
        root.destroy()
        return directory or None
    except Exception:
        return None

def _sync_from_text(state_key: str):
    st.session_state[state_key] = st.session_state.get(f"{state_key}_text", "").strip()

def dir_input(label: str, state_key: str, help_text: str | None = None, placeholder: str | None = None) -> str:
    """A text input with an adjacent 'Browse…' button that opens a folder picker."""
    if state_key not in st.session_state:
        st.session_state[state_key] = ""
    cols = st.columns([4, 1])
    with cols[0]:
        st.text_input(
            label,
            value=st.session_state[state_key],
            key=f"{state_key}_text",
            help=help_text,
            placeholder=placeholder,
            on_change=_sync_from_text,
            args=(state_key,),
        )
    with cols[1]:
        if st.button("Browse…", key=f"{state_key}_browse"):
            chosen = _pick_directory_dialog(initial=st.session_state[state_key] or None)
            if chosen:
                st.session_state[state_key] = chosen
                try:
                    st.rerun()
                except Exception:
                    st.experimental_rerun()
    return st.session_state[state_key]

# optional cache clear on start
if os.environ.get("CLEAR_CACHE_ON_START") == "1":
    try:
        st.cache_data.clear()
        st.cache_resource.clear()
    except Exception:
        pass

# ----------------------- Conserved regions (Merlin coordinates) -----------------------
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
    r'(?:_R?[12](?:_00[1-9])?)$', r'(?:_[12])$', r'(?:_sorted(?:_sorted)*)$',
    r'(?:_sort(?:ed)?)$', r'(?:_dedup(?:ed)?)$', r'(?:_mark(?:ed)?dups?)$',
    r'(?:_rmdup[s]?)$', r'(?:_trim(?:med)?(?:_paired)?)$', r'(?:_filtered)$',
    r'(?:\.bam)$',
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
    # type coercion for common numeric columns
    for col in ("avg_depth", "depth", "position", "total_reads", "total_hits"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out

@st.cache_data(show_spinner=False)
def load_all(csv_dir: str, blast_dir: str):
    gene_df     = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_gene.csv")))
    per_base_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_per_base.csv")))
    bam_df      = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_summary.csv")))
    blast_top5  = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_blast_top5.csv")))

    # BLAST raw hits (per-sample *_blast.tsv)
    blast_files = glob.glob(os.path.join(blast_dir, "*_blast.tsv"))
    hits = []
    for file in blast_files:
        try:
            df = pd.read_csv(file, sep="\t", header=None, names=[
                "qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore","stitle"
            ])
            if df.empty: continue
            df["sample"] = os.path.basename(file).split("_blast")[0]
            # numeric
            for c in ("pident","length","qstart","qend","sstart","send","evalue","bitscore"):
                df[c] = pd.to_numeric(df[c], errors="coerce")
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
    def _first_int(s):   m = re.search(r"(\d+)", s);  return int(m.group(1)) if m else 0
    def _first_float(s): m = re.search(r"(\d+(?:\.\d+)?)", s);  return float(m.group(1)) if m else 0.0

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
    flag_df = add_base_column(flag_df)
    return flag_df

@st.cache_data(show_spinner=False)
def load_vaf_tables(vaf_dir: str) -> pd.DataFrame:
    if not vaf_dir or not os.path.isdir(vaf_dir):
        return pd.DataFrame()
    paths = glob.glob(os.path.join(vaf_dir, "*_vaf.tsv"))
    rows = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t")
            if df.empty: continue
            sample = os.path.basename(p).replace("_vaf.tsv", "")
            df["sample"] = sample
            rows.append(df)
        except Exception as e:
            st.warning(f"Error reading VAF file {p}: {e}")
    if not rows:
        return pd.DataFrame()
    vaf = pd.concat(rows, ignore_index=True)
    for c in ("POS", "DP", "AD_ALT", "AF"):
        if c in vaf.columns:
            vaf[c] = pd.to_numeric(vaf[c], errors="coerce")
    vaf = add_base_column(vaf)
    return vaf

# ----------------------- Mixed strain & non-CMV heuristics -----------------------
def infer_mixed_strains(vaf_df: pd.DataFrame,
                        lo: float = 0.2, hi: float = 0.8,
                        count_likely: int = 10, frac_likely: float = 0.15,
                        count_possible: int = 3,  frac_possible: float = 0.05) -> str:
    """Heuristic: use mid-frequency variants (lo ≤ AF ≤ hi) to call mixed strains."""
    if vaf_df is None or vaf_df.empty or "AF" not in vaf_df.columns:
        return "Unknown"
    df = vaf_df.dropna(subset=["AF"])
    total = len(df)
    if total == 0:
        return "Unknown"
    mid = int(((df["AF"] >= lo) & (df["AF"] <= hi)).sum())
    frac = mid / total
    if mid >= count_likely and frac >= frac_likely:
        return "Likely"
    if mid >= count_possible and frac >= frac_possible:
        return "Possible"
    return "Unlikely"

def detect_non_cmv_species(blast_hits_all: pd.DataFrame, base: str) -> bool:
    """Detect notable non-CMV BLAST hits for a given base sample."""
    if blast_hits_all is None or blast_hits_all.empty or "base" not in blast_hits_all.columns:
        return False
    sub = blast_hits_all[blast_hits_all["base"] == base]
    if sub.empty:
        return False
    cmv_pat = re.compile(r"Human\s+herpesvirus\s*5|\bHCMV\b|Cytomegalovirus|NC_006273\.2", re.IGNORECASE)
    non = sub[~sub["stitle"].fillna("").str.contains(cmv_pat)]
    return len(non) >= 50  # adjustable threshold

# ----------------------- BLAST filters & E-value helpers -----------------------
def filter_blast_merlin(blast_hits: pd.DataFrame, min_identity=85.0,
                        min_length: int = 80, max_evalue: float = 1e-20) -> pd.DataFrame:
    """Keep CMV/Merlin hits and apply identity/length/E-value thresholds."""
    if blast_hits.empty:
        return blast_hits
    df = blast_hits.copy()
    cmv_mask = df["stitle"].str.contains(
        "Human herpesvirus 5|Human cytomegalovirus|Cytomegalovirus Merlin|NC_006273.2",
        case=False, na=False
    )
    df = df[cmv_mask]
    for c in ("pident","length","evalue","sstart","send"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[(df["pident"] >= float(min_identity)) &
            (df["length"] >= int(min_length)) &
            (df["evalue"] <= float(max_evalue))]
    return df

def _safe_evalues(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    x = np.where(~np.isfinite(x) | (x <= 0), 1e-300, x)  # treat zeros as tiny
    return x

def evalue_gmean(series: pd.Series) -> float:
    x = _safe_evalues(series)
    if x.size == 0:
        return np.nan
    return float(np.exp(np.mean(np.log(x))))

def evalue_neglog10(series: pd.Series) -> float:
    x = _safe_evalues(series)
    if x.size == 0:
        return np.nan
    return float(np.mean(-np.log10(x)))

def classify_significance(geom_e: float, strict: float, strong: float, weak: float) -> tuple[str,str]:
    """Return (label, comment) from geometric-mean E-value."""
    if not np.isfinite(geom_e):
        return ("No CMV hits", "No CMV-aligned reads in this set.")
    if geom_e <= strict:
        return ("Very unlikely due to chance", f"Geometric mean E={geom_e:.1e} ≤ {strict:.0e}.")
    if geom_e <= strong:
        return ("Unlikely due to chance", f"Geometric mean E={geom_e:.1e} ≤ {strong:.0e}.")
    if geom_e <= weak:
        return ("Weak/Indeterminate", f"Geometric mean E={geom_e:.1e} ≤ {weak:.0e}. Consider more data.")
    return ("Likely due to chance", f"Geometric mean E={geom_e:.1e} > {weak:.0e}.")

def map_hits_to_regions_merlin(df: pd.DataFrame) -> pd.DataFrame:
    """Duplicate hits for each conserved region they overlap (Merlin coords)."""
    if df is None or df.empty:
        return pd.DataFrame()
    merlin_mask = (df["sseqid"].astype(str).str.contains("NC_006273.2")) | (df["stitle"].astype(str).str.contains("Merlin", case=False))
    d = df[merlin_mask].copy()
    if d.empty:
        return pd.DataFrame()
    for c in ("sstart","send","evalue","pident","length"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    rows = []
    for _, r in d.dropna(subset=["sstart","send"]).iterrows():
        s0 = int(min(r["sstart"], r["send"]))
        s1 = int(max(r["sstart"], r["send"]))
        for reg in CONSERVE_REGIONS:
            if s1 >= reg["start"] and s0 <= reg["end"]:
                rr = dict(r)
                rr["region"] = reg["name"]
                rows.append(rr)
    return pd.DataFrame(rows)

def evalue_summary_tables(sel_hits_cmv: pd.DataFrame, bases: list[str],
                          th_strict: float, th_strong: float, th_weak: float):
    """Build (overall_df, per_region_df, messages_per_base)."""
    overall_rows, region_rows, messages = [], [], {}
    for b in bases:
        sub = sel_hits_cmv[sel_hits_cmv["base"] == b] if (sel_hits_cmv is not None and not sel_hits_cmv.empty) else pd.DataFrame()
        n = 0 if sub is None or sub.empty else int(len(sub))
        gmean_e = evalue_gmean(sub["evalue"]) if n > 0 else np.nan
        min_e   = float(np.min(_safe_evalues(sub["evalue"]))) if n > 0 else np.nan
        nl10    = evalue_neglog10(sub["evalue"]) if n > 0 else np.nan
        label, comment = classify_significance(gmean_e, th_strict, th_strong, th_weak)
        messages[b] = (label, comment, n)
        overall_rows.append({
            "base": b, "n_hits": n,
            "geo_mean_evalue": gmean_e, "min_evalue": min_e, "mean_-log10(E)": nl10,
            "significance": label, "comment": comment
        })
        sub_reg = map_hits_to_regions_merlin(sub)
        if not sub_reg.empty:
            for region, dfR in sub_reg.groupby("region"):
                r_n = int(len(dfR))
                r_gm = evalue_gmean(dfR["evalue"])
                r_min = float(np.min(_safe_evalues(dfR["evalue"])))
                r_nl10 = evalue_neglog10(dfR["evalue"])
                r_label, _ = classify_significance(r_gm, th_strict, th_strong, th_weak)
                region_rows.append({
                    "base": b, "region": region, "n_hits": r_n,
                    "geo_mean_evalue": r_gm, "min_evalue": r_min, "mean_-log10(E)": r_nl10,
                    "significance": r_label
                })
    return pd.DataFrame(overall_rows), pd.DataFrame(region_rows), messages

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
    # conserved region boxes
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
                      yaxis_title="Depth (log)", yaxis_type="log", xaxis_range=[0, 236000])
    return fig

# ----------------------- UI -----------------------
st.set_page_config(page_title="CMV Interactive Report", layout="wide")
st.title("📊 CMV Interactive Report")

with st.sidebar:
    st.header("Inputs")
    csv_dir = dir_input("CSV directory (gene/per_base/summary/blast_top5)",
                        state_key="csv_dir",
                        help_text="Folder with *_gene.csv, *_per_base.csv, *_summary.csv, *_blast_top5.csv")
    blast_dir = dir_input("BLAST directory (*_blast.tsv)",
                          state_key="blast_dir",
                          help_text="Folder with per-sample *_blast.tsv")
    vaf_dir = dir_input("VAF directory (*_vaf.tsv)", state_key="vaf_dir",
                        help_text="Folder with *_vaf.tsv (LoFreq)")
    bam_stats_dir = dir_input("BAM stats directory (*_flagstat.txt, *_stats.txt)",
                              state_key="bam_stats_dir",
                              help_text="Folder with samtools outputs")
    fastqc_dir = dir_input("FastQC directory (*_fastqc.html, *_fastqc.zip)",
                           state_key="fastqc_dir",
                           help_text="Folder with FastQC reports")

    # BLAST filtering knobs
    min_identity = st.number_input("BLAST min % identity (filter)", min_value=0.0, max_value=100.0, value=85.0, step=1.0)
    max_evalue   = st.number_input("BLAST max E-value (filter, ≤)", value=float(1e-20), format="%.1e")
    min_aln_len  = st.number_input("BLAST min alignment length (bp, filter)", min_value=0, value=80, step=10)
    min_cmv_reads = st.number_input("Presence call: min CMV hits (passing filters)", 0, 100000, 25, 1)

    with st.expander("Analysis thresholds (adjust as needed)", expanded=False):
        st.caption("E-value flags use geometric mean across hits.")
        # Alignment
        th_map_pct_min      = st.number_input("Min mapped %", 0.0, 100.0, 95.0, 1.0)
        th_proper_pct_min   = st.number_input("Min properly paired %", 0.0, 100.0, 90.0, 1.0)
        th_dup_pct_max      = st.number_input("Max duplicates %", 0.0, 100.0, 20.0, 1.0)
        th_singletons_max   = st.number_input("Max singletons %", 0.0, 100.0, 5.0, 0.5)
        # Coverage (no breadth)
        th_mean_depth_min   = st.number_input("Min mean depth per gene (×)", 0.0, 1e6, 30.0, 1.0)
        # Identity
        th_identity_min     = st.number_input("Min CMV % identity (mean)", 0.0, 100.0, 95.0, 0.5)
        # Mixed-strain heuristic
        th_mid_af_low       = st.number_input("Mixed AF lower bound", 0.0, 1.0, 0.2, 0.05)
        th_mid_af_high      = st.number_input("Mixed AF upper bound", 0.0, 1.0, 0.8, 0.05)
        th_mid_count_likely = st.number_input("Mixed: min mid-AF count (Likely)", 0, 100000, 10, 1)
        th_mid_frac_likely  = st.number_input("Mixed: min mid-AF fraction (Likely)", 0.0, 1.0, 0.15, 0.01)
        th_mid_count_poss   = st.number_input("Mixed: min mid-AF count (Possible)", 0, 100000, 3, 1)
        th_mid_frac_poss    = st.number_input("Mixed: min mid-AF fraction (Possible)", 0.0, 1.0, 0.05, 0.01)
        # E-value significance bands
        th_e_strict = st.number_input("E-value strict threshold (very unlikely)", value=float(1e-20), format="%.1e")
        th_e_strong = st.number_input("E-value strong threshold (unlikely)", value=float(1e-10), format="%.1e")
        th_e_weak   = st.number_input("E-value weak threshold (indeterminate)", value=float(1e-5),  format="%.1e")

# ----------------------- QC helpers -----------------------
def flag_alignment_sample(row, th_map_pct_min, th_proper_pct_min, th_dup_pct_max, th_singletons_max):
    issues = []
    if "mapped_pct" in row and pd.notna(row["mapped_pct"]) and row["mapped_pct"] < th_map_pct_min:
        issues.append(f"Low mapping rate ({row['mapped_pct']}% < {th_map_pct_min}%)")
    if "properly_paired_pct" in row and pd.notna(row["properly_paired_pct"]) and row["properly_paired_pct"] < th_proper_pct_min:
        issues.append(f"Low properly paired % ({row['properly_paired_pct']}% < {th_proper_pct_min}%)")
    if "duplicates_pct" in row and pd.notna(row["duplicates_pct"]) and row["duplicates_pct"] > th_dup_pct_max:
        issues.append(f"High duplicates ({row['duplicates_pct']}% > {th_dup_pct_max}%)")
    if "singletons_pct" in row and pd.notna(row["singletons_pct"]) and row["singletons_pct"] > th_singletons_max:
        issues.append(f"High singletons ({row['singletons_pct']}% > {th_singletons_max}%)")
    if "insert_size_mean" in row and pd.notna(row["insert_size_mean"]):
        mu = float(row["insert_size_mean"])
        if mu < 200 or mu > 600:
            issues.append(f"Unusual insert size mean (~{mu:.0f} bp)")
    return issues

def compute_gene_coverage_flags(gene_df, base, th_mean_depth_min):
    issues = []
    if gene_df is None or gene_df.empty:
        return issues
    g = gene_df[gene_df["base"] == base].copy()
    if g.empty: return issues
    low_depth_genes = g[(pd.to_numeric(g["avg_depth"], errors="coerce") < th_mean_depth_min)]
    if not low_depth_genes.empty:
        genes = ", ".join(sorted(low_depth_genes["gene"].astype(str).unique()))
        issues.append(f"Low mean depth (<{th_mean_depth_min}×) in gene(s): {genes}")
    return issues

def compute_identity_flags(blast_hits_df, base, th_identity_min):
    if blast_hits_df is None or blast_hits_df.empty:
        return [], None
    sub = blast_hits_df[blast_hits_df["base"] == base]
    if sub.empty or "pident" not in sub:
        return [], None
    nums = pd.to_numeric(sub["pident"], errors="coerce").dropna()
    mean_id = float(nums.mean()) if not nums.empty else None
    issues = []
    if mean_id is not None and mean_id < th_identity_min:
        issues.append(f"Mean CMV % identity low ({mean_id:.1f}% < {th_identity_min}%)")
    return issues, mean_id

def cmv_presence_call(sel_hits_cmv: pd.DataFrame, base: str, min_reads: int = 25):
    if sel_hits_cmv is None or sel_hits_cmv.empty or "base" not in sel_hits_cmv.columns:
        return "Not detected", 0
    sub = sel_hits_cmv[sel_hits_cmv["base"] == base]
    n = int(len(sub))
    if n >= int(min_reads): return "Detected", n
    if n > 0: return "Indeterminate", n
    return "Not detected", 0

# ----------------------- Main data loading & selection -----------------------
if csv_dir and blast_dir:
    gene_df, per_base_df, bam_df, blast_top5, blast_hits_all = load_all(csv_dir, blast_dir)
    align_df = load_alignment_stats(bam_stats_dir) if bam_stats_dir else pd.DataFrame()
    vaf_df   = load_vaf_tables(vaf_dir) if vaf_dir else pd.DataFrame()

    # Establish "core" bases from CMV CSVs; fall back to BLAST if needed
    core_sets = []
    for df in (gene_df, per_base_df, bam_df):
        if df is not None and not df.empty and "base" in df.columns:
            core_sets.append(set(df["base"].dropna().unique()))
    core_bases = set.union(*core_sets) if core_sets else set()
    if not core_bases:
        if blast_top5 is not None and not blast_top5.empty and "base" in blast_top5.columns:
            core_bases |= set(blast_top5["base"].dropna().unique())
        if blast_hits_all is not None and not blast_hits_all.empty and "base" in blast_hits_all.columns:
            core_bases |= set(blast_hits_all["base"].dropna().unique())
    all_bases = sorted(core_bases)

    def keep_bases(df):
        if df is None or df.empty or "base" not in df.columns: return df
        return df[df["base"].isin(core_bases)].copy()

    # Trim all inputs to "known" bases to keep everything consistent
    gene_df        = keep_bases(gene_df)
    per_base_df    = keep_bases(per_base_df)
    bam_df         = keep_bases(bam_df)
    blast_top5     = keep_bases(blast_top5)
    blast_hits_all = keep_bases(blast_hits_all)
    align_df       = keep_bases(align_df)
    vaf_df         = keep_bases(vaf_df)

    if not all_bases:
        st.warning("No samples found across inputs. Check your folder paths and file naming.")
        st.stop()

    with st.sidebar.expander("Debug: bases detected by source", expanded=False):
        def names(df): return sorted(df["base"].dropna().unique()) if (df is not None and not df.empty and "base" in df.columns) else []
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
            help="Pick one or more canonical sample IDs."
        )
        if not bases_selected:
            st.info("Select at least one base sample.")
            st.stop()

    # ---- Filter all tables to the selected bases (use f_* everywhere below) ----
    def _keep_selected(df):
        return df[df["base"].isin(bases_selected)].copy() if (df is not None and not df.empty and "base" in df.columns) else pd.DataFrame()

    f_align = _keep_selected(align_df)
    f_gene  = _keep_selected(gene_df)
    f_per   = _keep_selected(per_base_df)
    f_vaf   = _keep_selected(vaf_df)
    f_top5  = _keep_selected(blast_top5)

    # BLAST: restrict to CMV-ish titles and selected bases, then apply filters
    sel_hits_raw = (blast_hits_all[blast_hits_all["stitle"].str.contains(
        "Human herpesvirus 5|Human cytomegalovirus|Cytomegalovirus Merlin|NC_006273.2",
        case=False, na=False
    )].copy() if (blast_hits_all is not None and not blast_hits_all.empty) else pd.DataFrame())
    sel_hits_raw = sel_hits_raw[sel_hits_raw["base"].isin(bases_selected)] if not sel_hits_raw.empty else sel_hits_raw

    sel_hits = filter_blast_merlin(sel_hits_raw,
                                   min_identity=float(min_identity),
                                   min_length=int(min_aln_len),
                                   max_evalue=float(max_evalue))

    # ----------------------- 0) Analysis Summary -----------------------
    with st.expander("0. Analysis Summary", expanded=True):

        # helpers for summary
        def _agg_align_by_base(f_align: pd.DataFrame) -> pd.DataFrame:
            if f_align is None or f_align.empty:
                return pd.DataFrame()
            grp = f_align.groupby("base", as_index=False).agg({
                "total_reads":"sum",
                "mapped_reads":"sum",
                "properly_paired_reads":"sum",
                "paired_in_sequencing":"sum",
                "singletons":"sum",
                "duplicates":"sum",
                "insert_size_mean":"mean",
                "insert_size_sd":"mean"
            })
            grp["mapped_pct"] = np.where(grp["total_reads"]>0, grp["mapped_reads"]/grp["total_reads"]*100.0, np.nan)
            grp["properly_paired_pct"] = np.where(grp["paired_in_sequencing"]>0, grp["properly_paired_reads"]/grp["paired_in_sequencing"]*100.0, np.nan)
            grp["singletons_pct"] = np.where(grp["paired_in_sequencing"]>0, grp["singletons"]/grp["paired_in_sequencing"]*100.0, np.nan)
            grp["duplicates_pct"] = np.where(grp["total_reads"]>0, grp["duplicates"]/grp["total_reads"]*100.0, np.nan)
            return grp

        def _gene_depth_issues_by_base(f_gene: pd.DataFrame, th_mean_depth_min: float) -> dict:
            out = {}
            if f_gene is None or f_gene.empty:
                return out
            g_agg = f_gene.groupby(["base","gene"], as_index=False)["avg_depth"].mean(numeric_only=True)
            for b, gb in g_agg.groupby("base"):
                low = gb[pd.to_numeric(gb["avg_depth"], errors="coerce") < th_mean_depth_min]
                n = int(low.shape[0])
                names = list(low["gene"].astype(str).head(5))
                out[b] = (n, names)
            return out

        def _mean_identity_by_base(sel_hits_cmv: pd.DataFrame) -> pd.Series:
            if sel_hits_cmv is None or sel_hits_cmv.empty:
                return pd.Series(dtype=float)
            return sel_hits_cmv.groupby("base")["pident"].apply(lambda s: pd.to_numeric(s, errors="coerce").dropna().mean())

        rows = []
        align_agg = _agg_align_by_base(f_align)
        depth_issues = _gene_depth_issues_by_base(f_gene, th_mean_depth_min)
        mean_ident = _mean_identity_by_base(sel_hits)

        for b in bases_selected:
            arow = align_agg[align_agg["base"] == b].iloc[0] if (not align_agg.empty and (align_agg["base"]==b).any()) else None
            map_pct   = float(arow["mapped_pct"]) if arow is not None and pd.notna(arow["mapped_pct"]) else np.nan
            prop_pct  = float(arow["properly_paired_pct"]) if arow is not None and pd.notna(arow["properly_paired_pct"]) else np.nan
            dup_pct   = float(arow["duplicates_pct"]) if arow is not None and pd.notna(arow["duplicates_pct"]) else np.nan
            sing_pct  = float(arow["singletons_pct"]) if arow is not None and pd.notna(arow["singletons_pct"]) else np.nan

            align_fail_reasons = []
            if not np.isnan(map_pct)  and map_pct   < th_map_pct_min:    align_fail_reasons.append(f"mapped {map_pct:.1f}% < {th_map_pct_min}%")
            if not np.isnan(prop_pct) and prop_pct  < th_proper_pct_min: align_fail_reasons.append(f"proper {prop_pct:.1f}% < {th_proper_pct_min}%")
            if not np.isnan(dup_pct)  and dup_pct   > th_dup_pct_max:    align_fail_reasons.append(f"dup {dup_pct:.1f}% > {th_dup_pct_max}%")
            if not np.isnan(sing_pct) and sing_pct  > th_singletons_max: align_fail_reasons.append(f"singletons {sing_pct:.1f}% > {th_singletons_max}%")
            align_status = "OK" if not align_fail_reasons else "Issues: " + " · ".join(align_fail_reasons)

            n_low, names = depth_issues.get(b, (0, []))
            cov_status = "OK" if n_low == 0 else f"{n_low} gene(s) < {th_mean_depth_min}×" + (f" (e.g., {', '.join(names)})" if names else "")

            mid = mean_ident.get(b, np.nan)
            id_status = "OK" if (not np.isnan(mid) and mid >= th_identity_min) else ("No CMV hits" if np.isnan(mid) else f"Below threshold ({mid:.1f}% < {th_identity_min}%)")
            id_text = f"{mid:.1f}%" if not np.isnan(mid) else "n/a"

            sub = sel_hits[sel_hits["base"] == b] if not sel_hits.empty else pd.DataFrame()
            gmean_e = evalue_gmean(sub["evalue"]) if not sub.empty else np.nan
            e_label, _ = classify_significance(gmean_e, th_e_strict, th_e_strong, th_e_weak)
            e_col = f"{e_label} (geo-mean {gmean_e:.1e})" if np.isfinite(gmean_e) else "No CMV hits"

            non_cmv = detect_non_cmv_species(sel_hits_raw, b)
            non_cmv_status = "present" if non_cmv else "none"

            status, n_cmv = cmv_presence_call(sel_hits, b, min_reads=min_cmv_reads)
            presence_col = f"{status} (n={n_cmv})"

            # Mixed strains verdict for this base
            vf_b = f_vaf[f_vaf["base"] == b] if ("base" in f_vaf.columns and not f_vaf.empty) else pd.DataFrame()
            mixed_verdict = infer_mixed_strains(
                vf_b,
                lo=th_mid_af_low, hi=th_mid_af_high,
                count_likely=th_mid_count_likely, frac_likely=th_mid_frac_likely,
                count_possible=th_mid_count_poss,  frac_possible=th_mid_frac_poss
            )

            rows.append({
                "base": b,
                "CMV detected": presence_col,
                "Alignment": align_status if align_status == "OK" else f"⚠️ {align_status}",
                "Coverage (mean per gene)": "OK" if cov_status == "OK" else f"⚠️ {cov_status}",
                "CMV % identity (mean)": f"{id_status} ({id_text})" if id_status != "OK" else f"OK ({id_text})",
                "E-value significance (overall)": e_col,
                "Mixed strains": mixed_verdict,
                "Non-CMV hits": non_cmv_status,
            })

        summary_df = pd.DataFrame(rows)
        st.dataframe(summary_df, use_container_width=True)

    # ----------------------- 1) Alignment -----------------------
    with st.expander("1. Alignment Statistics (BWA-MEM)", expanded=True):
        st.subheader("QC flags")
        if not f_align.empty:
            for _, row in f_align.iterrows():
                issues = flag_alignment_sample(row, th_map_pct_min, th_proper_pct_min, th_dup_pct_max, th_singletons_max)
                sample = row.get("sample", "unknown")
                if issues:
                    st.warning(f"🧬 **{sample}** — " + " · ".join(issues))
                else:
                    st.success(f"🧬 **{sample}** — Within thresholds")

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
            st.info("No alignment stats found for selected samples.")

    # ----------------------- 2. Average Depth per Gene -----------------------
    with st.expander("2. Average Depth per Gene", expanded=True):
        st.subheader("QC flags")
        if not f_gene.empty:
            for b in bases_selected:
                issues = compute_gene_coverage_flags(f_gene, base=b, th_mean_depth_min=th_mean_depth_min)
                if issues: st.warning(f"📉 **{b}** — " + " · ".join(issues))
                else:      st.success(f"📉 **{b}** — Mean depth per gene within thresholds")

            st.dataframe(f_gene.drop(columns=["base"]), use_container_width=True)
            fig = plot_avg_depth(f_gene)
            if fig: st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No gene-level coverage data for selected samples.")

    # ----------------------- 3. BLAST Summary, % Identity & E-values -----------------------
    with st.expander("3. BLAST Summary, % Identity & E-value Significance", expanded=True):
        st.subheader("QC flags — % identity")
        if not sel_hits.empty:
            for b in bases_selected:
                issues, mean_id = compute_identity_flags(sel_hits, base=b, th_identity_min=th_identity_min)
                if issues:
                    more = f" (mean={mean_id:.1f}%)" if mean_id is not None else ""
                    st.warning(f"🧪 **{b}** — " + " · ".join(issues) + more)
                else:
                    if mean_id is not None:
                        st.success(f"🧪 **{b}** — Mean CMV identity OK ({mean_id:.1f}%)")
        else:
            st.info("No CMV-filtered BLAST hits after current filters for selected samples.")

        st.subheader("E-value significance (overall & by conserved region)")
        overall_tbl, region_tbl, messages = evalue_summary_tables(
            sel_hits_cmv=sel_hits, bases=bases_selected,
            th_strict=th_e_strict, th_strong=th_e_strong, th_weak=th_e_weak
        )

        # Friendly banners (no inline ternary that leaks DeltaGenerator)
        for b in bases_selected:
            label, comment, n = messages.get(b, ("No CMV hits", "No CMV-aligned reads.", 0))
            prefix = f"🔬 **{b}** — "
            msg = f"{prefix}{label}: {comment} (n={n} hits)"
            if label.startswith("Very"):
                st.success(msg)
            elif label.startswith("Unlikely"):
                st.info(msg)
            elif label.startswith("Weak"):
                st.warning(msg)
            elif label.startswith("Likely"):
                st.warning(msg)
            else:  # "No CMV hits"
                st.info(msg)

        if not overall_tbl.empty:
            disp = overall_tbl.copy()
            disp["geo_mean_evalue"] = disp["geo_mean_evalue"].apply(lambda x: f"{x:.1e}" if pd.notna(x) else "n/a")
            disp["min_evalue"]      = disp["min_evalue"].apply(lambda x: f"{x:.1e}" if pd.notna(x) else "n/a")
            disp["mean_-log10(E)"]  = disp["mean_-log10(E)"].apply(lambda x: f"{x:.2f}" if pd.notna(x) else "n/a")
            st.markdown("**Overall E-value summary (selected samples)**")
            st.dataframe(disp, use_container_width=True)

        if not region_tbl.empty:
            dispR = region_tbl.copy().sort_values(["base","region"])
            dispR["geo_mean_evalue"] = dispR["geo_mean_evalue"].apply(lambda x: f"{x:.1e}" if pd.notna(x) else "n/a")
            dispR["min_evalue"]      = dispR["min_evalue"].apply(lambda x: f"{x:.1e}" if pd.notna(x) else "n/a")
            dispR["mean_-log10(E)"]  = dispR["mean_-log10(E)"].apply(lambda x: f"{x:.2f}" if pd.notna(x) else "n/a")
            st.markdown("**Per conserved region E-value summary (Merlin coordinates; selected samples)**")
            st.dataframe(dispR, use_container_width=True)
        else:
            st.info("No Merlin-specific region hits available for current filters.")

        if not f_top5.empty:
            st.subheader("Top BLAST species per sample (selected)")
            st.dataframe(f_top5.drop(columns=["base"]), use_container_width=True)
            if "species" in f_top5.columns and "total_hits" in f_top5.columns:
                fig = px.bar(f_top5, x="sample", y="total_hits", color="species",
                             title="Top BLAST Species per Sample and Read")
                fig.update_xaxes(tickangle=45)
                st.plotly_chart(fig, use_container_width=True)

        if not sel_hits.empty:
            st.subheader("BLAST % identity across genome (selected, filtered)")
            fig = plot_blast_full_heatmap(sel_hits, bin_size=100)
            if fig: st.plotly_chart(fig, use_container_width=True)

            st.subheader("BLAST % identity by conserved region (selected, filtered)")
            tabs = st.tabs([r["name"] for r in CONSERVE_REGIONS])
            for tab, region in zip(tabs, CONSERVE_REGIONS):
                with tab:
                    fig = plot_blast_region_heatmap(sel_hits, region=region)
                    if fig:
                        st.plotly_chart(fig, use_container_width=True)

    # ----------------------- 4. Per-base Coverage Heatmaps -----------------------
    with st.expander("4. Per-base Coverage Heatmaps", expanded=False):
        if not f_per.empty:
            tabs = st.tabs([r["name"] for r in CONSERVE_REGIONS])
            for tab, region in zip(tabs, CONSERVE_REGIONS):
                with tab:
                    fig = plot_per_base_heatmap(f_per, gene=region["name"], bin_size=None)
                    if fig: st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No per-base depth data for selected samples.")

    # ----------------------- 5. Genome-wide Coverage -----------------------
    with st.expander("5. Genome-wide Coverage", expanded=False):
        if not f_per.empty:
            fig = plot_genome_coverage(f_per)
            if fig: st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No per-base depth data for genome-wide view (selected).")

    # ----------------------- 6. FastQC Reports -----------------------
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
                try:
                    sample = os.path.basename(z).replace("_fastqc.zip", "")
                    with zipfile.ZipFile(z) as zf:
                        inner = [n for n in zf.namelist() if n.endswith("/summary.txt")]
                        if inner:
                            with zf.open(inner[0]) as fh:
                                for line in io.TextIOWrapper(fh, encoding="utf-8", errors="ignore"):
                                    parts = line.strip().split("\t")
                                    if len(parts) >= 3:
                                        status, module, _filename = parts[0], parts[1], parts[2]
                                        all_sum.append({"sample": sample, "module": module, "status": status})
                except Exception as e:
                    st.warning(f"Failed to read {z}: {e}")
            if all_sum:
                df_sum = pd.DataFrame(all_sum)
                df_sum["base"] = df_sum["sample"].apply(canonical_sample)
                df_sum = df_sum[df_sum["base"].isin(bases_selected)].drop(columns=["base"])
                st.subheader("FastQC module summary (selected)")
                pivot = df_sum.pivot_table(index=["sample"], columns="module", values="status",
                                           aggfunc=lambda x: ", ".join(sorted(set(x))))
                st.dataframe(pivot.fillna("N/A"), use_container_width=True)

            reports = html_files + zip_files
            if reports:
                choice = st.selectbox("Preview a FastQC report", options=["(select)"] + reports)
                if choice and choice != "(select)":
                    st.markdown(f"**Preview:** {os.path.basename(choice)}")
                    try:
                        if choice.endswith(".zip"):
                            with zipfile.ZipFile(choice) as zf:
                                html_files = [f for f in zf.namelist() if f.endswith("fastqc_report.html")]
                                if html_files:
                                    with zf.open(html_files[0]) as f:
                                        html_content = f.read().decode("utf-8", errors="ignore")
                                        components.html(html_content, height=800, scrolling=True)
                        else:
                            with open(choice, "r", encoding="utf-8", errors="ignore") as f:
                                html_content = f.read()
                                components.html(html_content, height=800, scrolling=True)
                    except Exception as e:
                        st.warning(f"Preview failed for {choice}: {e}")
            else:
                st.info("No FastQC HTML or ZIP files found for selected samples.")
        else:
            st.info("Provide a valid FastQC directory to list reports.")

    # ----------------------- 7. Variants (LoFreq VAF) -----------------------
    with st.expander("7. Variants (LoFreq VAF)", expanded=True):
        if f_vaf is None or f_vaf.empty:
            st.info("No VAF tables loaded for selected samples.")
        else:
            st.subheader("QC flags — mixed strains")
            for b in bases_selected:
                vf_b = f_vaf[f_vaf["base"] == b] if ("base" in f_vaf.columns) else pd.DataFrame()
                verdict = infer_mixed_strains(
                    vf_b,
                    lo=th_mid_af_low, hi=th_mid_af_high,
                    count_likely=th_mid_count_likely, frac_likely=th_mid_frac_likely,
                    count_possible=th_mid_count_poss,  frac_possible=th_mid_frac_poss
                )
                msg = f"🧬 **{b}** — Mixed HCMV strains: **{verdict}**"
                if verdict == "Likely":
                    st.warning(msg)
                elif verdict == "Possible":
                    st.info(msg)
                else:
                    st.success(msg)

            st.subheader("VAF table (selected samples)")
            st.dataframe(
                f_vaf[[c for c in ["sample","CHROM","POS","REF","ALT","DP","AD_ALT","AF","base"] if c in f_vaf.columns]]
                .sort_values(["sample","POS"]).reset_index(drop=True)
                .drop(columns=["base"]),
                use_container_width=True,
            )

            c1, c2 = st.columns(2)
            with c1:
                fig = go.Figure()
                for sample, sdf in f_vaf.groupby("sample"):
                    sdf = sdf.sort_values("POS")
                    fig.add_trace(go.Scatter(x=sdf["POS"], y=sdf["AF"], mode="markers",
                                             name=sample, hovertemplate="POS=%{x}<br>AF=%{y:.3f}<extra>"))
                fig.update_layout(title="Variant Allele Frequency (AF) across genome",
                                  xaxis_title="Position (bp)", yaxis_title="AF (0-1)",
                                  yaxis=dict(range=[0,1]))
                st.plotly_chart(fig, use_container_width=True)
            with c2:
                fig = px.histogram(f_vaf, x="AF", color="sample", nbins=40, title="Distribution of AF values")
                fig.update_xaxes(range=[0,1])
                st.plotly_chart(fig, use_container_width=True)

else:
    st.info("Enter at least CSV and BLAST directories to load data.")
