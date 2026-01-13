"""
CMV Interactive Report (Streamlit)

This script builds a Streamlit web app for reviewing Cytomegalovirus (CMV)
sequencing pipeline outputs across one or more samples.

It is designed to pull together several outputs produced by a diagnostic
sequencing workflow, including:

- Gene-level depth summaries (*_gene.csv)
- Per-base depth tables (*_per_base.csv) (optionally loaded, because these can
  be large)
- Alignment summary metrics (e.g., cmv_alignment_summary.csv or *_summary.csv)
- BLAST summaries (e.g., *_blast_top5.csv) and raw BLAST hits (*_blast.tsv)
- LoFreq variant allele frequency tables (*_vaf.tsv)
- samtools alignment QC outputs (*_flagstat.txt, *_stats.txt)
- FastQC reports (*_fastqc.html, *_fastqc.zip)
- snpEff-annotated VCFs (*.annot.vcf / *.annot.vcf.gz) to flag potential
  antiviral resistance mutations in UL97 and UL54
"""

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

# ----------------------- Folder picker (Browse…) -----------------------
def _pick_directory_dialog(initial: str | None = None) -> str | None:
    """
    Open a native folder-selection dialog and return the chosen directory path.

    This function uses `tkinter` to open the operating system's folder picker. 
    It is mainly used to support a "Browse…" button in the Streamlit app,
    so the user can choose a folder rather than typing a path.

    Args:
        initial (str | None):
            The directory path to open the dialog in initially.
            - If a valid path is provided, the dialog starts there.
            - If None, the dialog defaults to the current working directory
              (`os.getcwd()`).

    Returns:
        str | None:
            The selected directory path as a string if the user chooses one.
            Returns None if:
            - the user cancels the folder picker,
            - the dialog fails to open,
            - or tkinter is unavailable / an exception occurs.
    """
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        try:
            root.call("wm", "attributes", ".", "-topmost", True)
        except Exception:
            pass
        directory = filedialog.askdirectory(initialdir=initial or os.getcwd())
        root.destroy()
        return directory or None
    except Exception:
        return None



def dir_input(
    label: str,
    state_key: str,
    help_text: str | None = None,
    placeholder: str | None = None,
) -> str:
    """
    Create a Streamlit directory input made of:
    1) a text box to type a path, and
    2) a "Browse…" button to pick a folder using a native dialog.

    Flow:
    - If a pending value exists, apply it to st.session_state[state_key].
    - Show a text input for manual path entry.
    - Show a Browse button:
        - If clicked, open a folder picker (tkinter).
        - If a folder is chosen, store it under the pending key and call
          st.rerun() so the text input updates cleanly.
    - Return the final path string (with whitespace stripped).

    Args:
        label (str):
            Label displayed next to/above the text input in the Streamlit UI.
            Example: "CSV directory".

        state_key (str):
            The Streamlit session_state key used to store the text input value.
            This must be unique per widget in the app (e.g., "csv_dir").

        help_text (str | None):
            Optional help tooltip shown in the UI (the little ⓘ tooltip).

        placeholder (str | None):
            Optional placeholder string shown inside the text input when empty.

    Returns:
        str:
            The directory path currently stored for this widget, with leading
            and trailing whitespace removed.
            Returns an empty string ("") if nothing has been entered/selected.
    """
    pending_key = f"{state_key}_pending"
    if pending_key in st.session_state:
        st.session_state[state_key] = st.session_state.pop(pending_key)

    cols = st.columns([4, 1])
    with cols[0]:
        st.text_input(
            label,
            key=state_key,
            help=help_text,
            placeholder=placeholder,
        )
    with cols[1]:
        if st.button("Browse…", key=f"{state_key}_browse"):
            chosen = _pick_directory_dialog(initial=st.session_state.get(state_key) or os.getcwd())
            if chosen:
                st.session_state[pending_key] = chosen
                st.rerun()

    return st.session_state.get(state_key, "").strip()



# ----------------------- Cache behaviour -----------------------
# Optional: clear Streamlit caches when the app starts/reruns.
# This is mainly for debugging, so you can force the app to reload data fresh
# instead of reusing cached results from st.cache_data / st.cache_resource.
# It only runs if the environment variable CLEAR_CACHE_ON_START is set to "1".
if os.environ.get("CLEAR_CACHE_ON_START") == "1":
    try:
        # Clear cached data results (e.g. DataFrames returned from @st.cache_data functions).
        st.cache_data.clear()
        # Clear cached resources (e.g. expensive objects stored with @st.cache_resource).
        st.cache_resource.clear()
    except Exception:
        # If cache clearing is not supported or fails for any reason, do nothing
        # so the app still runs.
        pass


# ----------------------- Common read options -----------------------
_COMMON_READ_KW = dict(engine="c", low_memory=False)

# ----------------------- Conserved regions (Merlin coordinates) -----------------------
# Genomic coordinates (Merlin reference) for regions of interest.
# These are used later to:
# - label/highlight important CMV genes on plots (heatmaps, genome coverage)
# - define which positions to extract for region-specific visualisations
CONSERVE_REGIONS = [
    {"name": "IE1", "start": 172328, "end": 174090},
    {"name": "UL54", "start": 78193, "end": 81922},
    {"name": "gpUL55", "start": 82065, "end": 84789},
    {"name": "UL75", "start": 109223, "end": 111452},
    {"name": "UL97", "start": 141797, "end": 143921},
    {"name": "UL83", "start": 120655, "end": 122341},
]

# ----------------------- Standardising sample names -----------------------
# Patterns to strip from sample names to get standardised base IDs.
_SUFFIX_PATTERNS = [
    r"(?:_R?[12](?:_00[1-9])?)$",
    r"(?:_[12])$",
    r"(?:_sorted(?:_sorted)*)$",
    r"(?:_sort(?:ed)?)$",
    r"(?:_dedup(?:ed)?)$",
    r"(?:_mark(?:ed)?dups?)$",
    r"(?:_rmdup[s]?)$",
    r"(?:_trim(?:med)?(?:_paired)?)$",
    r"(?:_filtered)$",
    r"(?:\.bam)$",
]


def canonical_sample(name: str) -> str:
    """
    Standardise a sample identifier by removing common pipeline/file suffixes.

    The function repeatedly applies each regular expression in `_SUFFIX_PATTERNS`
    until no more changes occur. 

    Examples:
        - "SAMPLE001_R1" -> "SAMPLE001"
        - "SAMPLE001_sorted_sorted" -> "SAMPLE001"
        - "SAMPLE001_R2_dedup.bam" -> "SAMPLE001"

    Args:
        name (str):
            Input sample name or filename-derived sample ID.

    Returns:
        str:
            Canonical/base sample ID with known suffixes removed. Any trailing
            "_" or "-" characters left behind after stripping are also removed.

    Notes:
        - Matching is case-insensitive (re.IGNORECASE).
        - The suffix patterns are defined in the module-level `_SUFFIX_PATTERNS`
          list.
    """
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
    """
    Add a 'base' column to the DataFrame by standardising the 'sample' column.

    Args:
        df (pd.DataFrame):
            Input DataFrame. Expected to contain a column named 'sample' where
            each value is a sample ID (often derived from a filename).

    Returns:
        pd.DataFrame:
            A DataFrame with an additional 'base' column if possible.
            - If 'sample' exists and the DataFrame is not empty, a copy of the
              DataFrame is returned with the new 'base' column added.
            - Otherwise, the original DataFrame is returned unchanged.

    """
    if df is not None and not df.empty and "sample" in df.columns:
        df = df.copy()
        df["base"] = df["sample"].astype(str).apply(canonical_sample)
    return df


# ----------------------- Data loaders -----------------------

def _read_csv_list(files, **read_kwargs) -> pd.DataFrame:
    """
    Read and concatenate multiple CSV/TSV files into a single DataFrame.
    
    This is used throughout the app to load pipeline outputs that are
    produced per-sample (e.g. "*_gene.csv", "*_summary.csv", "*_blast_top5.csv").
    It loops over a list of file paths, reads each one with pandas, and then
    concatenates them into a single table.

    Args:
        files (Iterable[str]):
            List of file paths to read. Each file is expected to be in CSV
            format (or TSV if specified in `read_kwargs`).

        **read_kwargs:
            Additional keyword arguments to pass to `pd.read

    Returns:
        pd.DataFrame:
            A single DataFrame containing the concatenated data from all
            successfully read files.

    Notes:
        - Files that are empty (0 bytes) are skipped to avoid warnings and wasted IO.
        - If a file exists but contains no parsable data, it is also skipped.
    """
    dfs = []
    for f in files:
        try:
            # Skip empty (0-byte) files to avoid warnings and wasted IO.
            if os.path.getsize(f) == 0:
                continue

            # Read file with shared defaults + any caller-specific parameters.
            df = pd.read_csv(f, **_COMMON_READ_KW, **read_kwargs)
            if not df.empty:
                dfs.append(df)

        except pd.errors.EmptyDataError:
            # File exists but contains no rows/columns that pandas can parse.
            continue

        except Exception as e:
            # Keep the app running even if one file fails to read.
            st.warning(f"Error reading {f}: {e}")

    # If nothing was loaded successfully, return an empty DataFrame.
    if not dfs:
        return pd.DataFrame()

    # Combine all loaded DataFrames into one.
    out = pd.concat(dfs, ignore_index=True)

    # Type coercion for common numeric columns (vectorized).
    # This prevents string/object dtype from breaking numeric thresholds/plots.
    for col in ("avg_depth", "depth", "position", "total_reads", "total_hits"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    # CMV alignment metric columns (if present).
    for col in (
        "cmv_mapped_reads_mapq",
        "cmv_rpm_mapq",
        "cmv_breadth_mapq",
        "cmv_cov_bases_mapq",
        "cmv_ref_bases",
        "total_reads_in_cmv_bam",
    ):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    return out



@st.cache_data(ttl=3600, show_spinner=False)
def load_all(csv_dir: str, blast_dir: str | None):
    """
    Load the main result tables used by the Streamlit CMV report.

    This function gathers the pipeline outputs from the provided input
    directories and returns them as pandas DataFrames. It is cached using
    Streamlit's `@st.cache_data` to avoid repeatedly re-reading files on each
    rerun of the app (Streamlit reruns the script frequently when widgets
    change).

    What gets loaded:
        1) gene_df: Concatenated gene-level coverage tables from "*_gene.csv".
        2) per_base_df (lightweight placeholder):Instead of loading full "*_per_base.csv" content 
            (which can be very large), this builds a small DataFrame containing only sample names
            that have a per-base file present. 
        3) bam_df:Alignment summary table. If a combined file "cmv_alignment_summary.csv" exists, 
            it is used. Otherwise, it falls back to concatenating "*_summary.csv" files.
        4) blast_top5: Concatenated BLAST summary tables from "*_blast_top5.csv".
        5) blast_hits: Optional - concatenated raw BLAST hits from "*_blast.tsv" (outfmt 6).
            The app keeps identity/length/coordinate/bitscore/stitle columns 

    After loading, a "base" column is added to tables that contain a "sample" column.
    This allows the app to group multiple related filenames (e.g., R1/R2, sorted/deduped outputs) 
    under one base sample ID.

    Args:
        csv_dir (str):
            Directory containing the main CSV outputs from the pipeline.
        blast_dir (str | None):
            Optional directory containing per-sample BLAST hit tables

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            A 5-tuple in this order:
              - gene_df: gene-level coverage data
              - per_base_df: placeholder DataFrame listing samples with per-base files
              - bam_df: alignment summary data
              - blast_top5: top BLAST species summary data
              - blast_hits: raw BLAST hits (E-value removed), may be empty
    """
    # ------------------- 1) Gene-level coverage -------------------
    gene_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_gene.csv")))

    # ------------------- 2) Per-base depth (placeholder only) -------------------
    # This avoids loading large per-base tables unless the user enables those sections.
    per_base_files = glob.glob(os.path.join(csv_dir, "*_per_base.csv"))
    per_base_df = pd.DataFrame(
        {"sample": [os.path.basename(p).replace("_per_base.csv", "") for p in per_base_files]}
    )
    per_base_df = add_base_column(per_base_df)

    # ------------------- 3) Alignment summary -------------------
    # Prefer one combined summary file if present; otherwise use per-sample summaries.
    cmv_summary = os.path.join(csv_dir, "cmv_alignment_summary.csv")
    if os.path.isfile(cmv_summary):
        bam_df = _read_csv_list([cmv_summary])
    else:
        bam_df = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_summary.csv")))

    # ------------------- 4) BLAST top species summary -------------------
    blast_top5 = _read_csv_list(glob.glob(os.path.join(csv_dir, "*_blast_top5.csv")))

    # ------------------- 5) BLAST raw hits (optional) -------------------
    # BLAST raw hits (per-sample *_blast.tsv). We keep pident/length/coords/bitscore/stitle.
    blast_hits = pd.DataFrame()
    if blast_dir and os.path.isdir(blast_dir):
        blast_files = glob.glob(os.path.join(blast_dir, "*_blast.tsv"))
        hits = []
        if blast_files:
            for file in blast_files:
                try:
                    # Skip 0-byte BLAST files to avoid warning spam and wasted IO.
                    if os.path.getsize(file) == 0:
                        continue

                    # BLAST outfmt 6 commonly includes evalue and stitle.
                    df = pd.read_csv(
                        file,
                        sep="\t",
                        header=None,
                        names=[
                            "qseqid",
                            "sseqid",
                            "pident",
                            "length",
                            "mismatch",
                            "gapopen",
                            "qstart",
                            "qend",
                            "sstart",
                            "send",
                            "evalue",  # read but NOT used; dropped below
                            "bitscore",
                            "stitle",
                        ],
                        **_COMMON_READ_KW,
                    )
                    if df.empty:
                        continue

                    # Infer sample name from filename prefix.
                    df["sample"] = os.path.basename(file).split("_blast")[0]

                    # Coerce numeric columns so thresholds/plots work reliably.
                    for c in (
                        "pident",
                        "length",
                        "mismatch",
                        "gapopen",
                        "qstart",
                        "qend",
                        "sstart",
                        "send",
                        "bitscore",
                    ):
                        if c in df.columns:
                            df[c] = pd.to_numeric(df[c], errors="coerce")

                    # Keep ID/text columns as strings for consistency.
                    if "qseqid" in df.columns:
                        df["qseqid"] = df["qseqid"].astype("string")
                    if "sseqid" in df.columns:
                        df["sseqid"] = df["sseqid"].astype("string")
                    if "stitle" in df.columns:
                        df["stitle"] = df["stitle"].astype("string")

                    if "evalue" in df.columns:
                        df = df.drop(columns=["evalue"])

                    hits.append(df)

                except pd.errors.EmptyDataError:
                    # File exists but has no parseable content.
                    continue
                except Exception as e:
                    # Non-fatal: warn in the UI but continue processing.
                    st.warning(f"Error reading {file}: {e}")

            blast_hits = pd.concat(hits, ignore_index=True) if hits else pd.DataFrame()

    # ------------------- Attach standardised base IDs -------------------
    gene_df = add_base_column(gene_df)
    bam_df = add_base_column(bam_df)
    blast_top5 = add_base_column(blast_top5)
    blast_hits = add_base_column(blast_hits)

    return gene_df, per_base_df, bam_df, blast_top5, blast_hits


@st.cache_data(ttl=3600, show_spinner=False)
def load_per_base_for_selected(csv_dir: str, bases: tuple[str, ...]) -> pd.DataFrame:
    """
    Load per-base depth tables (*_per_base.csv) for the selected base sample IDs.

    The app only reads them when the user enables sections that require per-base
    plots (e.g., coverage heatmaps or genome-wide depth plots).

    For each base ID, the function searches for files matching:
        "{base}*_per_base.csv"

    The function then:
    - de-duplicates the file list, sorts it, reads and concatenates them,
    - adds the standard "base" column for consistent grouping/filtering.

    Args:
        csv_dir (str):
            Directory containing per-base depth CSV files (pattern: "*_per_base.csv").

        bases (tuple[str, ...]):
            Tuple of canonical base sample IDs selected in the UI.
            Example: ("SAMPLE001", "SAMPLE002")

    Returns:
        pd.DataFrame:
            Combined per-base depth DataFrame for the requested bases.
            Typically contains columns such as:
              - sample (original sample ID from filenames or file content)
              - gene
              - position
              - depth
              - base (added by add_base_column)

            Returns an empty DataFrame if no matching files are found or all
            matched files are empty/unreadable.
    """
    files = []
    for b in bases:
        # Match per-base files that start with the base ID.
        files.extend(glob.glob(os.path.join(csv_dir, f"{b}*_per_base.csv")))
    df = _read_csv_list(sorted(set(files)))
    df = add_base_column(df)
    return df



@st.cache_data(ttl=3600, show_spinner=False)
def load_alignment_stats(bam_stats_dir: str) -> pd.DataFrame:
    """
    Extract samtools alignment QC outputs into a single table.

    This function reads:
      - "*_flagstat.txt" (samtools flagstat): totals, mapped %, duplicates,
        properly paired %, singletons, etc.
      - "*_stats.txt" (samtools stats): insert size mean/SD and read length mean
        from the "SN" summary lines.

    It returns one row per sample and adds a canonical "base" column so results
    can be grouped by base sample ID in the app.

    Args:
        bam_stats_dir (str):
            Directory containing samtools outputs:
            "*_flagstat.txt" and "*_stats.txt".

    Returns:
        pd.DataFrame:
            Combined alignment metrics per sample. 

            Returns an empty DataFrame if no files are found/parsed successfully.
    """
    def _first_int(s):
        m = re.search(r"(\d+)", s)
        return int(m.group(1)) if m else 0

    def _first_float(s):
        m = re.search(r"(\d+(?:\.\d+)?)", s)
        return float(m.group(1)) if m else 0.0

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
            flag_rows.append(
                {
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
            )
        except Exception as e:
            st.warning(f"Failed to parse flagstat {f}: {e}")

    for s in glob.glob(os.path.join(bam_stats_dir, "*_stats.txt")):
        try:
            sample = os.path.basename(s).replace("_stats.txt", "")
            ins_mean = ins_sd = read_len = None
            with open(s, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith("SN"):
                        continue
                    if "insert size average:" in line:
                        ins_mean = _first_float(line)
                    elif "insert size standard deviation:" in line:
                        ins_sd = _first_float(line)
                    elif "average length:" in line:
                        read_len = _first_float(line)
            sn_rows.append(
                {"sample": sample, "insert_size_mean": ins_mean, "insert_size_sd": ins_sd, "read_length_mean": read_len}
            )
        except Exception as e:
            st.warning(f"Failed to parse stats {s}: {e}")

    flag_df = pd.DataFrame(flag_rows)
    sn_df = pd.DataFrame(sn_rows)
    if not flag_df.empty and not sn_df.empty:
        flag_df = flag_df.merge(sn_df, on="sample", how="left")
    if "sample" in flag_df.columns:
        flag_df = flag_df.sort_values("sample").reset_index(drop=True)
    flag_df = add_base_column(flag_df)
    return flag_df


@st.cache_data(ttl=3600, show_spinner=False)
def load_vaf_tables(vaf_dir: str) -> pd.DataFrame:
    """
    Load and combine LoFreq VAF tables from a directory.

    This function looks for files ending in "*_vaf.tsv", reads them into pandas,
    adds a "sample" column based on the filename, and concatenates everything
    into one DataFrame. 


    Args:
        vaf_dir (str):
            Path to a directory containing LoFreq-style VAF tables named
            "*_vaf.tsv".

    Returns:
        pd.DataFrame:
            Combined VAF table across all samples. Typical columns include:
            CHROM, POS, REF, ALT, DP, AD_ALT, AF, sample, and derived base.
            Returns an empty DataFrame if nothing is loaded.
    """
    if not vaf_dir or not os.path.isdir(vaf_dir):
        return pd.DataFrame()
    paths = glob.glob(os.path.join(vaf_dir, "*_vaf.tsv"))
    rows = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t", **_COMMON_READ_KW)
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
    for c in ("POS", "DP", "AD_ALT", "AF"):
        if c in vaf.columns:
            vaf[c] = pd.to_numeric(vaf[c], errors="coerce")
    vaf = add_base_column(vaf)
    return vaf



# ----------------------- Antiviral resistance parsing (snpEff VCF) -----------------------
UL97_COMMON = {"A594V", "L595S", "M460I", "M460V", "C592G", "H520Q", "C603W"}
UL97_HOTSPOTS = set([460, 520] + list(range(590, 608)))  # inclusive 590–607

UL54_NAMED = {"V812L", "A834P"}  # exact named calls
UL54_HOTSPOT_CODONS = {812, 834}  # if other nonsyn changes occur here, flag as hotspot
UL54_DEL981_982_PAT = re.compile(r"(?:^|[^A-Za-z])(?:del|Del|DEL)\s*9?81[_\-–]?\s*982(?:[^A-Za-z]|$)|Δ\s*981\s*[_\-–]?\s*982")
HGVS_P_PAT = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|Ter)")
ONE_LETTER = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
    "Ter": "*",
}


def _hgvs_p_to_short(hgvs_p: str) -> tuple[str | None, int | None, str | None]:
    """
    Convert a snpEff HGVS protein string into a short (one-letter) representation.

    Example:
        "p.Ala594Val" -> ("A", 594, "V")

    The function uses a regular expression to extract
      - 3-letter amino acid code (from)
      - position (integer)
      - 3-letter amino acid code (to), "=" (synonymous), or "Ter" (stop)

    If parsing fails, it returns (None, None, None).

    Args:
        hgvs_p (str):
            HGVS protein annotation string (e.g., "p.Ala594Val").

    Returns:
        tuple[str | None, int | None, str | None]:
            (aa_from, position, aa_to) where:
              - aa_from is a one-letter AA code or None
              - position is an int or None
              - aa_to is a one-letter AA code (or "*" for stop) or None
    """
    if not hgvs_p:
        return (None, None, None)
    m = HGVS_P_PAT.search(hgvs_p)
    if not m:
        return (None, None, None)
    aa_from_3, pos_str, aa_to_3 = m.group(1), m.group(2), m.group(3)
    pos = int(pos_str)
    f = ONE_LETTER.get(aa_from_3, None)
    t = ONE_LETTER.get(aa_to_3, None)
    return (f, pos, t)



def _ann_entries(info: str) -> list[list[str]]:
    """
    Extract and split snpEff ANN entries from a VCF INFO field.

    snpEff writes annotations into the INFO field under the key "ANN=".
    This function:
      1) locates the ANN=... part of INFO
      2) splits multiple annotations by comma
      3) splits each annotation into fields using '|'

    The returned structure is a list of lists, where each inner list is the
    pipe-split ANN fields.

    Args:
        info (str):
            The INFO column string from a VCF record.

    Returns:
        list[list[str]]:
            A list of ANN entries split on '|'.
            Returns an empty list if ANN is not present or INFO is empty.
    """
    if not info:
        return []
    ann_raw = None
    for kv in info.split(";"):
        if kv.startswith("ANN="):
            ann_raw = kv[4:]
            break
    if not ann_raw:
        return []
    entries = ann_raw.split(",")
    return [e.split("|") for e in entries if e]



@st.cache_data(ttl=3600, show_spinner=False)
def load_resistance_flags_from_vcfs(vcf_dir_annot: str) -> pd.DataFrame:
    """
    Read snpEff-annotated VCF files and summarise antiviral resistance findings.

    This function scans "*.annot.vcf" / "*.annot.vcf.gz" files and looks for
    UL97 and UL54 protein changes in ANN annotations. Based on a small set of
    rules (common named mutations + hotspot logic + a specific deletion pattern),
    it assigns per-base resistance flags and generates a human-readable summary.

    Output is one row per base sample with:
      - resistance_flag: True if any qualifying finding exists
      - details: string summary (multiple findings joined by " | ")

    Args:
        vcf_dir_annot (str):
            Directory containing snpEff annotated VCFs:
            "*.annot.vcf" or "*.annot.vcf.gz".

    Returns:
        pd.DataFrame:
            Columns: ["base", "resistance_flag", "details"].
            Empty (with those columns) if directory/files are missing or no
            qualifying findings are detected.
    """
    if not vcf_dir_annot or not os.path.isdir(vcf_dir_annot):
        return pd.DataFrame(columns=["base", "resistance_flag", "details"])
    vcfs = sorted(glob.glob(os.path.join(vcf_dir_annot, "*.annot.vcf*")))
    if not vcfs:
        return pd.DataFrame(columns=["base", "resistance_flag", "details"])

    findings = {}

    def _open_any(path):
        """
        Open a text handle for a VCF, supporting both plain and gzipped files.

        Args:
            path (str):
                File path to a VCF (.vcf or .vcf.gz).

        Returns:
            TextIO:
                Open file handle in text mode.
        """
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    for vcf in vcfs:
        sample_guess = os.path.basename(vcf)
        sample_guess = re.sub(r"\.annot\.vcf(?:\.gz)?$", "", sample_guess)
        base = canonical_sample(sample_guess)
        try:
            with _open_any(vcf) as fh:
                for line in fh:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    info = parts[7]
                    for ann in _ann_entries(info):
                        gene = ann[3] if len(ann) > 3 else ""
                        hgvs_p = ann[10] if len(ann) > 10 else ""
                        gene = (gene or "").strip().upper()
                        hgvs_p = (hgvs_p or "").strip()
                        if not gene or not hgvs_p:
                            continue

                        aa_from, pos, aa_to = _hgvs_p_to_short(hgvs_p)

                        # UL97 rules
                        if gene == "UL97":
                            label = None
                            if aa_from and aa_to and pos:
                                short = f"{aa_from}{pos}{aa_to}"
                                if short in UL97_COMMON:
                                    label = f"UL97 GCV/VGCV resistance ({short})"
                                elif pos in UL97_HOTSPOTS and aa_to != "*" and aa_to != aa_from:
                                    label = f"UL97 hotspot {pos} ({short}; GCV/VGCV resistance likely)"
                            if label:
                                findings.setdefault(base, set()).add(label)

                        # UL54 rules
                        elif gene == "UL54":
                            if aa_from and pos:
                                short = f"{aa_from}{pos}{aa_to}" if aa_to else None

                                if short and short in UL54_NAMED:
                                    findings.setdefault(base, set()).add(f"UL54 {short} (GCV/CDV/PFA cross-resistance)")

                                if (hgvs_p and (("981" in hgvs_p and "982" in hgvs_p) and (("del" in hgvs_p.lower()) or ("Δ" in hgvs_p)))) \
                                   or (UL54_DEL981_982_PAT.search(hgvs_p or "")):
                                    findings.setdefault(base, set()).add("UL54 Δ981–982 (GCV+CDV+PFA resistance)")

                                if aa_to and aa_to != "*" and aa_to != aa_from and pos in UL54_HOTSPOT_CODONS:
                                    findings.setdefault(base, set()).add(f"UL54 hotspot {pos} variant ({short}); verify resistance")

        except Exception as e:
            st.warning(f"Failed to parse {vcf}: {e}")

    rows = []
    for b, flags in findings.items():
        rows.append({"base": b, "resistance_flag": True, "details": " | ".join(sorted(flags))})
    if not rows:
        return pd.DataFrame(columns=["base", "resistance_flag", "details"])
    df = pd.DataFrame(rows).sort_values("base").reset_index(drop=True)
    return df



# ----------------------- Mixed strain & non-CMV heuristics -----------------------
def infer_mixed_strains(
    vaf_df: pd.DataFrame,
    lo: float = 0.2,
    hi: float = 0.8,
    count_likely: int = 10,
    frac_likely: float = 0.15,
    count_possible: int = 3,
    frac_possible: float = 0.05,
) -> str:
    """
    Infer likelihood of mixed CMV strains from VAF distribution.
    Rationale:
      - Single-strain infections tend to have variants near 0 or near 1.
      - Mixed strains tend to produce more variants with allele frequencies in
        the middle range (e.g., 0.2–0.8).

    This function counts variants where lo <= AF <= hi and then returns:
      - "Likely" if both count and fraction thresholds are met
      - "Possible" if smaller thresholds are met
      - "Unlikely" otherwise
      - "Unknown" if input is missing/empty/AF not available

    Args:
        vaf_df (pd.DataFrame):
            VAF table for a single base sample (or already-filtered subset).
            Must contain an "AF" column (allele frequency in 0–1).
        lo (float):
            Lower bound for "mid" AF.
        hi (float):
            Upper bound for "mid" AF.
        count_likely (int):
            Minimum number of mid-AF variants for a "Likely" call.
        frac_likely (float):
            Minimum fraction of mid-AF variants for a "Likely" call.
        count_possible (int):
            Minimum number of mid-AF variants for a "Possible" call.
        frac_possible (float):
            Minimum fraction of mid-AF variants for a "Possible" call.

    Returns:
        str:
            One of: "Likely", "Possible", "Unlikely", or "Unknown".
    """
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
    """
    Flag whether a sample has a notable number of non-CMV BLAST hits.

    The function:
      - selects BLAST hits for a given base sample
      - removes hits whose stitle looks CMV-like (HCMV / HHV-5 / Merlin)
      - returns True if the remaining (non-CMV) hits count is >= 10

    Args:
        blast_hits_all (pd.DataFrame):
            Combined BLAST hits table (typically created in load_all()).
            Expected columns include "base" and "stitle".
        base (str):
            Standard base sample ID to check.

    Returns:
        bool:
            True if >= 10 non-CMV hits are found, otherwise False.
    """
    if blast_hits_all is None or blast_hits_all.empty or "base" not in blast_hits_all.columns:
        return False
    sub = blast_hits_all[blast_hits_all["base"] == base]
    if sub.empty:
        return False
    cmv_pat = re.compile(r"Human\s+herpesvirus\s*5|\bHCMV\b|Cytomegalovirus|NC_006273\.2", re.IGNORECASE)
    non = sub[~sub["stitle"].fillna("").str.contains(cmv_pat)]
    return len(non) >= 10  # adjustable threshold



# ----------------------- BLAST filters (NO E-VALUE) -----------------------
def filter_blast_merlin(blast_hits: pd.DataFrame, min_identity=85.0, min_length: int = 80) -> pd.DataFrame:
    """
    Filter BLAST hits to CMV/Merlin-like titles and basic quality thresholds.

    This function keeps only hits whose "stitle" suggests CMV then applies:
      - pident >= min_identity
      - length >= min_length


    Args:
        blast_hits (pd.DataFrame):
            Raw/combined BLAST hits. Expected columns include:
            "stitle", "pident", "length", "sstart", "send".
        min_identity (float):
            Minimum percent identity to keep (0–100).
        min_length (int):
            Minimum alignment length to keep (bp).

    Returns:
        pd.DataFrame:
            Filtered BLAST hits table. Returns an empty DataFrame if input is
            empty or None.
    """
    if blast_hits is None or blast_hits.empty:
        return pd.DataFrame()

    df = blast_hits.copy()
    cmv_mask = df["stitle"].astype(str).str.contains(
        "Human herpesvirus 5|Human cytomegalovirus|Cytomegalovirus Merlin|NC_006273.2",
        case=False,
        na=False,
    )
    df = df[cmv_mask]

    for c in ("pident", "length", "sstart", "send"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df[(df["pident"] >= float(min_identity)) & (df["length"] >= int(min_length))]
    return df



# ----------------------- Plotting -----------------------
def plot_alignment_stacked(df: pd.DataFrame):
    """
    Create a stacked horizontal bar chart of mapped vs unmapped reads per sample.

    Args:
        df (pd.DataFrame):
            Alignment stats table with columns: "sample", "mapped_reads",
            "unmapped_reads".

    Returns:
        plotly.graph_objects.Figure | None:
            Plotly Figure if data exists, otherwise None.
    """
    if df.empty:
        return None
    fig = go.Figure()
    fig.add_bar(y=df["sample"], x=df["mapped_reads"], orientation="h", name="Mapped reads")
    fig.add_bar(y=df["sample"], x=df["unmapped_reads"], orientation="h", name="Unmapped reads")
    fig.update_layout(barmode="stack", xaxis_title="Reads", yaxis_title="Sample", title="Mapped vs Unmapped Reads per Sample")
    return fig


def plot_alignment_pct(df: pd.DataFrame, col: str, title: str):
    """
    Create a bar chart for a percentage-based alignment metric.

    Args:
        df (pd.DataFrame):
            Alignment stats table with a "sample" column and the requested metric column.
        col (str):
            Column name to plot on the y-axis (e.g., "mapped_pct").
        title (str):
            Plot title.

    Returns:
        plotly.graph_objects.Figure | None:
            Plotly Figure if column exists and data is non-empty, otherwise None.
    """
    if df.empty or col not in df:
        return None
    fig = px.bar(df, x="sample", y=col, title=title)
    fig.update_yaxes(title="%")
    fig.update_xaxes(tickangle=45)
    return fig


def plot_avg_depth(gene_df: pd.DataFrame):
    """
    Create a bar chart of average depth per CMV gene, coloured by sample.

    Args:
        gene_df (pd.DataFrame):
            Gene-level depth table with columns: "gene", "avg_depth", "sample".

    Returns:
        plotly.graph_objects.Figure | None:
            Plotly Figure if data exists, otherwise None.
    """
    if gene_df.empty:
        return None
    fig = px.bar(gene_df, x="gene", y="avg_depth", color="sample", title="Average Depth per CMV Gene")
    fig.update_yaxes(title="Average Depth")
    fig.update_xaxes(title="Gene")
    return fig


def plot_blast_full_heatmap(blast_hits: pd.DataFrame, bin_size=100):
    """
    Build a genome-wide BLAST % identity heatmap in fixed-size bins.

    Each BLAST hit contributes its % identity to bins across the hit's subject
    coordinate range (sstart..send). Bins are aggregated by mean % identity.

    Conserved regions are overlaid as outlined rectangles.

    Args:
        blast_hits (pd.DataFrame):
            Filtered BLAST hits with columns: "sample", "sstart", "send", "pident".
        bin_size (int):
            Bin size in bp for genome aggregation.

    Returns:
        plotly.graph_objects.Figure | None:
            Heatmap Figure, or None if inputs do not produce any binned rows.
    """
    if blast_hits.empty:
        return None
    rows = []
    for _, r in blast_hits.iterrows():
        s0 = int(min(r["sstart"], r["send"]))
        s1 = int(max(r["sstart"], r["send"]))
        b0 = s0 - (s0 % bin_size)
        for b in range(b0, max(b0 + bin_size, s1), bin_size):
            rows.append([r["sample"], b, float(r["pident"])])
    if not rows:
        return None
    heat_df = pd.DataFrame(rows, columns=["sample", "position", "pident"])
    pivot = heat_df.groupby(["sample", "position"]).pident.mean().unstack(fill_value=0).sort_index(axis=1)
    fig = px.imshow(
        pivot,
        aspect="auto",
        color_continuous_scale="RdBu",
        zmin=0,
        zmax=100,
        labels=dict(color="% Identity"),
        title=f"BLAST % Identity Across Genome ({bin_size} bp bins)",
    )
    shapes = []
    for reg in CONSERVE_REGIONS:
        x0 = reg["start"] - (reg["start"] % bin_size)
        x1 = reg["end"] - (reg["end"] % bin_size)
        shapes.append(dict(type="rect", xref="x", yref="paper", x0=x0, x1=x1, y0=0, y1=1, line=dict(color="black"), fillcolor="rgba(0,0,0,0)"))
    fig.update_layout(shapes=shapes)
    return fig


def plot_blast_region_heatmap(blast_hits: pd.DataFrame, region: dict):
    """
    Build a per-base-position BLAST % identity heatmap for a single region.

    For each sample, the function averages % identity across all hits covering
    each position in the region (start..end).

    Args:
        blast_hits (pd.DataFrame):
            Filtered BLAST hits with columns: "sample", "sstart", "send", "pident".
        region (dict):
            Region dictionary with keys: "name", "start", "end".

    Returns:
        plotly.graph_objects.Figure | None:
            Heatmap Figure for the region, or None if inputs are empty.
    """
    if blast_hits.empty:
        return None
    rstart, rend = int(region["start"]), int(region["end"])
    region_len = rend - rstart + 1
    columns = np.arange(rstart, rend + 1)
    samples = sorted(blast_hits["sample"].unique())
    sums = np.zeros((len(samples), region_len), dtype=float)
    counts = np.zeros_like(sums, dtype=int)
    sidx = {s: i for i, s in enumerate(samples)}
    for _, r in blast_hits.iterrows():
        s0 = int(min(r["sstart"], r["send"]))
        s1 = int(max(r["sstart"], r["send"]))
        left, right = max(s0, rstart), min(s1, rend)
        if left > right:
            continue
        i = sidx.get(r["sample"])
        a, b = left - rstart, right - rstart + 1
        p = float(r["pident"])
        sums[i, a:b] += p
        counts[i, a:b] += 1
    means = np.divide(sums, counts, out=np.zeros_like(sums, dtype=float), where=counts > 0)
    pivot = pd.DataFrame(means, index=samples, columns=columns)
    fig = px.imshow(
        pivot,
        aspect="auto",
        color_continuous_scale="RdBu",
        zmin=0,
        zmax=100,
        labels=dict(color="% Identity"),
        title=f"BLAST % Identity: {region['name']} ({rstart}-{rend})",
    )
    return fig


def plot_per_base_heatmap(per_base_df: pd.DataFrame, gene: str, bin_size: int | None = None):
    """
    Plot per-base depth as a heatmap for a specific gene/region.

    The region coordinates are taken from CONSERVE_REGIONS. If bin_size is
    provided, positions are grouped into bins and mean depth is plotted.

    Args:
        per_base_df (pd.DataFrame):
            Per-base depth table with columns including: "sample", "gene",
            "position", "depth".
        gene (str):
            Region name (must match a CONSERVE_REGIONS "name" entry).
        bin_size (int | None):
            If provided and > 1, depth is averaged within bins of this size.

    Returns:
        plotly.graph_objects.Figure | None:
            Heatmap Figure, or None if there is no data for the region.
    """
    region_map = {r["name"]: (r["start"], r["end"]) for r in CONSERVE_REGIONS}
    start, end = region_map[gene]
    data = per_base_df[(per_base_df["gene"] == gene) & (per_base_df["position"] >= start) & (per_base_df["position"] <= end)].copy()
    if data.empty:
        return None
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
    """
    Plot genome-wide per-base depth for each sample as line traces (log y-axis).

    Conserved regions are marked with shaded vertical rectangles.

    Args:
        per_base_df (pd.DataFrame):
            Per-base depth table with columns: "sample", "position", "depth".

    Returns:
        plotly.graph_objects.Figure | None:
            Line plot Figure, or None if input is empty.
    """
    if per_base_df.empty:
        return None
    fig = go.Figure()
    for sample, df in per_base_df.groupby("sample"):
        df = df.sort_values("position")
        fig.add_trace(go.Scatter(x=df["position"], y=df["depth"], mode="lines", name=sample))
    for region in CONSERVE_REGIONS:
        fig.add_vrect(
            x0=region["start"],
            x1=region["end"],
            fillcolor="orange",
            opacity=0.2,
            line_width=0,
            annotation_text=region["name"],
            annotation_position="top",
        )
    fig.update_layout(
        title="Per-base Depth Across Genome",
        xaxis_title="Genomic Position",
        yaxis_title="Depth (log)",
        yaxis_type="log",
        xaxis_range=[0, 236000],
    )
    return fig

# ----------------------- UI -----------------------
# Configure the Streamlit page (tab title + wide layout for plots/tables).
st.set_page_config(page_title="CMV Interactive Report", layout="wide")

# Main page title shown at the top of the app.
st.title("📊 CMV Interactive Report")

# Everything inside this block appears in the left-hand sidebar.
with st.sidebar:
    st.header("Inputs")

    # Folder containing the main pipeline CSV outputs (gene/per-base/summary/BLAST top5).
    csv_dir = dir_input(
        "CSV directory (gene/per_base/summary/blast_top5)",
        state_key="csv_dir",
        help_text="Folder with *_gene.csv, *_per_base.csv, *_summary.csv, *_blast_top5.csv",
    )

    # Optional folder containing raw BLAST hit tables for identity plots.
    blast_dir = dir_input(
        "BLAST directory (*_blast.tsv) [optional]",
        state_key="blast_dir",
        help_text="Folder with per-sample *_blast.tsv (optional; for BLAST identity plots)",
    )

    # Folder containing LoFreq variant tables (VAF).
    vaf_dir = dir_input(
        "VAF directory (*_vaf.tsv)",
        state_key="vaf_dir",
        help_text="Folder with *_vaf.tsv (LoFreq)",
    )

    # Folder containing samtools QC outputs (flagstat + stats).
    bam_stats_dir = dir_input(
        "BAM stats directory (*_flagstat.txt, *_stats.txt)",
        state_key="bam_stats_dir",
        help_text="Folder with samtools outputs",
    )

    # Folder containing FastQC reports (either HTML files or ZIP archives).
    fastqc_dir = dir_input(
        "FastQC directory (*_fastqc.html, *_fastqc.zip)",
        state_key="fastqc_dir",
        help_text="Folder with FastQC reports",
    )

    # Folder containing snpEff-annotated VCFs used for resistance flagging.
    vcf_annot_dir = dir_input(
        "Annotated VCF directory (*.annot.vcf[.gz])",
        state_key="vcf_annot_dir",
        help_text="Output of snpEff step (from annotate_vcfs_snpeff.sh)",
    )

    # BLAST filtering settings (identity + alignment length; no E-values).
    min_identity = st.number_input(
        "BLAST min % identity (filter)",
        min_value=0.0,
        max_value=100.0,
        value=85.0,
        step=1.0,
    )
    min_aln_len = st.number_input(
        "BLAST min alignment length (bp, filter)",
        min_value=0,
        value=80,
        step=10,
    )

    # Advanced thresholds (collapsed by default to keep the sidebar tidy).
    with st.expander("Analysis thresholds (adjust as needed)", expanded=False):
        # Alignment QC thresholds.
        th_map_pct_min = st.number_input("Min mapped %", 0.0, 100.0, 95.0, 1.0)
        th_proper_pct_min = st.number_input("Min properly paired %", 0.0, 100.0, 90.0, 1.0)
        th_dup_pct_max = st.number_input("Max duplicates %", 0.0, 100.0, 20.0, 1.0)
        th_singletons_max = st.number_input("Max singletons %", 0.0, 100.0, 5.0, 0.5)

        # Coverage threshold for gene-level depth.
        th_mean_depth_min = st.number_input("Min mean depth per gene (×)", 0.0, 1e6, 30.0, 1.0)

        # Identity threshold for mean CMV % identity (from BLAST).
        th_identity_min = st.number_input("Min CMV % identity (mean)", 0.0, 100.0, 95.0, 0.5)

        # Mixed strain heuristic thresholds based on mid-range AF variants.
        th_mid_af_low = st.number_input("Mixed AF lower bound", 0.0, 1.0, 0.2, 0.05)
        th_mid_af_high = st.number_input("Mixed AF upper bound", 0.0, 1.0, 0.8, 0.05)
        th_mid_count_likely = st.number_input("Mixed: min mid-AF count (Likely)", 0, 100000, 10, 1)
        th_mid_frac_likely = st.number_input("Mixed: min mid-AF fraction (Likely)", 0.0, 1.0, 0.15, 0.01)
        th_mid_count_poss = st.number_input("Mixed: min mid-AF count (Possible)", 0, 100000, 3, 1)
        th_mid_frac_poss = st.number_input("Mixed: min mid-AF fraction (Possible)", 0.0, 1.0, 0.05, 0.01)

        # CMV detection thresholds based on alignment HQ mapped reads.
        st.markdown("**CMV detection (alignment HQ mapped)**")
        th_detect_mapped = st.number_input("Detected threshold (mapped_hq ≥)", 0.0, 1e12, 57.0, 1.0)
        th_indeterminate_low = st.number_input("Indeterminate lower bound (≥)", value=float(1e-12), format="%.1e")
        th_not_detected_high = st.number_input("Not detected upper bound (<)", value=float(1e-12), format="%.1e")

        # Thresholds for whether UL97/UL54 have enough coverage to assess resistance.
        st.markdown("**Antiviral resistance assessability (UL97/UL54)**")
        th_res_gene_depth_min = st.number_input(
            "Min mean depth in UL97/UL54 to call 'Not detected' (×)",
            0.0,
            1e6,
            30.0,
            1.0,
        )
        th_res_gene_min_count = st.number_input("Min number of assessable genes (UL97/UL54)", 0, 2, 1, 1)

    # Section toggles so the user can disable slower plots/tables.
    with st.expander("Sections to show", expanded=True):
        st.caption("Turn off heavier sections to speed things up.")

        # High-level sections.
        show_alignment = st.checkbox("1. Alignment Statistics (BWA-MEM)", value=True)
        show_depth_per_gene = st.checkbox("2. Average Depth per Gene", value=True)
        show_blast = st.checkbox("3. BLAST / % Identity (no E-values)", value=True)

        # Only show BLAST sub-options if BLAST section is enabled.
        if show_blast:
            st.markdown("**3.x BLAST sub-sections**")
            show_blast_bar = st.checkbox("3.a Top BLAST species bar chart", value=True)
            show_blast_heatmap = st.checkbox("3.b Genome % identity heatmap", value=False)
            show_blast_regions = st.checkbox("3.c Per conserved region heatmaps", value=False)
        else:
            # Keep values defined even if BLAST is turned off.
            show_blast_bar = show_blast_heatmap = show_blast_regions = False

        # Heavier per-base plots are off by default to reduce load time.
        show_per_base_heatmaps = st.checkbox("4. Per-base Coverage Heatmaps", value=False)
        show_genome_wide = st.checkbox("5. Genome-wide Coverage", value=False)

        # FastQC preview and variant plots.
        show_fastqc = st.checkbox("6. FastQC Reports", value=True)
        show_variants = st.checkbox("7. Variants (LoFreq VAF)", value=True)

    # Data loading can be expensive, so it is gated behind a button.
    # The app only reads files once the user clicks "Load / Refresh".
    load_refresh_clicked = st.button("Load / Refresh", key="load_refresh_button")


# Build a simple "signature" of the current input directory paths.
# If any directory changes, we reset the load gate so the user must re-load.
_dirs_sig = (csv_dir, blast_dir, vaf_dir, bam_stats_dir, fastqc_dir, vcf_annot_dir)

# If the signature has changed since the last run, mark data as not committed/loaded.
if st.session_state.get("_dirs_sig") != _dirs_sig:
    st.session_state["_dirs_sig"] = _dirs_sig
    st.session_state["_load_committed"] = False

# If the user clicks the button, commit the load so downstream code can run.
if load_refresh_clicked:
    st.session_state["_load_committed"] = True



# ----------------------- QC helpers (used by summary table only) -----------------------
def flag_alignment_sample(row, th_map_pct_min, th_proper_pct_min, th_dup_pct_max, th_singletons_max):
    """
    Generate human-readable QC issue strings for a single alignment row.

    Args:
        row (pd.Series | dict-like):
            A row containing alignment metrics (mapped_pct, duplicates_pct, etc.).
        th_map_pct_min (float):
            Minimum acceptable mapped percentage.
        th_proper_pct_min (float):
            Minimum acceptable properly paired percentage.
        th_dup_pct_max (float):
            Maximum acceptable duplicates percentage.
        th_singletons_max (float):
            Maximum acceptable singletons percentage.

    Returns:
        list[str]:
            List of QC issue messages. Empty list means no issues detected.
    """
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
    """
    Identify low-depth genes for a given base sample.

    Args:
        gene_df (pd.DataFrame):
            Gene-level coverage table with columns including "base", "gene", "avg_depth".
        base (str):
            Canonical base sample ID to evaluate.
        th_mean_depth_min (float):
            Minimum acceptable mean depth (×).

    Returns:
        list[str]:
            Issue messages describing genes below the depth threshold.
    """
    issues = []
    if gene_df is None or gene_df.empty:
        return issues
    g = gene_df[gene_df["base"] == base].copy()
    if g.empty:
        return issues
    low_depth_genes = g[(pd.to_numeric(g["avg_depth"], errors="coerce") < th_mean_depth_min)]
    if not low_depth_genes.empty:
        genes = ", ".join(sorted(low_depth_genes["gene"].astype(str).unique()))
        issues.append(f"Low mean depth (<{th_mean_depth_min}×) in gene(s): {genes}")
    return issues


def compute_identity_flags(blast_hits_df, base, th_identity_min):
    """
    Compute mean CMV % identity for a base sample and flag if below threshold.

    Args:
        blast_hits_df (pd.DataFrame):
            BLAST hits table containing "base" and "pident".
        base (str):
            Canonical base sample ID to evaluate.
        th_identity_min (float):
            Minimum acceptable mean % identity.

    Returns:
        tuple[list[str], float | None]:
            (issues, mean_identity)
            - issues: list of messages (empty if OK/no data)
            - mean_identity: float mean identity if calculable, else None
    """
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


def _agg_cmv_mapped_hq_by_base(bam_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate CMV high-quality mapped reads per base sample.

    This sums "cmv_mapped_reads_mapq" across rows grouped by "base" and returns
    a table with "base" and "cmv_mapped_hq".

    Args:
        bam_df (pd.DataFrame):
            Alignment summary table expected to contain "base" and
            "cmv_mapped_reads_mapq".

    Returns:
        pd.DataFrame:
            Columns: ["base", "cmv_mapped_hq"].
            Returns an empty DataFrame with these columns if inputs are missing.
    """
    if bam_df is None or bam_df.empty or "base" not in bam_df.columns:
        return pd.DataFrame(columns=["base", "cmv_mapped_hq"])
    if "cmv_mapped_reads_mapq" not in bam_df.columns:
        return pd.DataFrame(columns=["base", "cmv_mapped_hq"])
    grp = bam_df.groupby("base", as_index=False).agg({"cmv_mapped_reads_mapq": "sum"})
    grp = grp.rename(columns={"cmv_mapped_reads_mapq": "cmv_mapped_hq"})
    grp["cmv_mapped_hq"] = pd.to_numeric(grp["cmv_mapped_hq"], errors="coerce").fillna(0.0)
    return grp


def cmv_presence_call_mapped_hq(
    cmv_mapped_hq: float,
    detected_threshold: float = 57.0,
    indeterminate_low: float = 1e-12,
    not_detected_high: float = 1e-12,
) -> str:
    """
    Call CMV presence status based on high-quality mapped reads.

    Logic:
      - "Not detected"   if x < not_detected_high
      - "Indeterminate"  if indeterminate_low <= x < detected_threshold
      - "Detected"       if x >= detected_threshold

    Args:
        cmv_mapped_hq (float):
            High-quality mapped reads (or a summed value) for CMV.
        detected_threshold (float):
            Threshold at or above which CMV is called "Detected".
        indeterminate_low (float):
            Lower bound for the "Indeterminate" band.
        not_detected_high (float):
            Upper bound for "Not detected" (exclusive).

    Returns:
        str:
            One of: "Not detected", "Indeterminate", "Detected".
    """
    x = float(cmv_mapped_hq) if cmv_mapped_hq is not None and np.isfinite(cmv_mapped_hq) else 0.0
    if x < float(not_detected_high):
        return "Not detected"
    if x < float(detected_threshold) and x >= float(indeterminate_low):
        return "Indeterminate"
    return "Detected"


def resistance_assessability_by_base(gene_df: pd.DataFrame, base: str, depth_min: float) -> tuple[bool, dict]:
    """
    Decide whether UL97 and/or UL54 have enough coverage to assess resistance.

    This does not call resistance. It only checks whether mean depth in UL97/UL54
    is >= depth_min (using gene-level avg_depth). The function reports how many
    of the two genes are "assessable".

    Args:
        gene_df (pd.DataFrame):
            Gene-level depth table containing "base", "gene", "avg_depth".
        base (str):
            Canonical base sample ID to evaluate.
        depth_min (float):
            Minimum mean depth to consider a gene assessable.

    Returns:
        tuple[bool, dict]:
            (assessable, detail)
            - assessable: True if at least one of UL97/UL54 meets depth_min
            - detail: dict with keys:
                "UL97_mean_depth", "UL54_mean_depth", "genes_assessable"
    """
    detail = {"UL97_mean_depth": np.nan, "UL54_mean_depth": np.nan, "genes_assessable": 0}
    if gene_df is None or gene_df.empty or "base" not in gene_df.columns:
        return False, detail
    sub = gene_df[gene_df["base"] == base]
    if sub.empty or "gene" not in sub.columns or "avg_depth" not in sub.columns:
        return False, detail

    g = sub.copy()
    g["gene_u"] = g["gene"].astype(str).str.upper()
    g["avg_depth"] = pd.to_numeric(g["avg_depth"], errors="coerce")

    ul97 = g[g["gene_u"] == "UL97"]["avg_depth"].dropna()
    ul54 = g[g["gene_u"] == "UL54"]["avg_depth"].dropna()

    ul97_mean = float(ul97.mean()) if not ul97.empty else np.nan
    ul54_mean = float(ul54.mean()) if not ul54.empty else np.nan
    detail["UL97_mean_depth"] = ul97_mean
    detail["UL54_mean_depth"] = ul54_mean

    genes_ok = 0
    if np.isfinite(ul97_mean) and ul97_mean >= depth_min:
        genes_ok += 1
    if np.isfinite(ul54_mean) and ul54_mean >= depth_min:
        genes_ok += 1
    detail["genes_assessable"] = genes_ok
    return genes_ok > 0, detail

# ----------------------- Main data loading & selection -----------------------
# This block runs only after:
#   1) the user has provided a CSV directory (csv_dir), AND
#   2) the user has clicked the "Load / Refresh" button (stored as _load_committed)
#
# This is the “load everything, figure out which samples exist,
# let the user pick samples, then prepare filtered dataframes for the plots” step.
if csv_dir and st.session_state.get("_load_committed", False):
    # Load the main outputs produced by the pipeline.
    # Each of these is typically a pandas DataFrame, one row per sample ("base") or per metric.
    gene_df, per_base_df, bam_df, blast_top5, blast_hits_all = load_all(csv_dir, blast_dir)

    # Optional inputs (only loaded if the corresponding directory path was provided).
    align_df = load_alignment_stats(bam_stats_dir) if bam_stats_dir else pd.DataFrame()
    vaf_df = load_vaf_tables(vaf_dir) if vaf_dir else pd.DataFrame()
    resist_df = load_resistance_flags_from_vcfs(vcf_annot_dir) if vcf_annot_dir else pd.DataFrame()

    # Establish "core" bases from CMV CSVs; fall back to BLAST if needed
    #
    # "base" is your canonical sample identifier used across outputs.
    # Here we build a set of all sample IDs we can find across the loaded tables.
    core_sets = []
    for df in (gene_df, per_base_df, bam_df):
        # Only attempt to pull sample IDs if:
        # - df exists and has rows
        # - it contains a "base" column (your canonical sample key)
        if df is not None and not df.empty and "base" in df.columns:
            core_sets.append(set(df["base"].dropna().unique()))

    # Merge all sets of sample IDs into one.
    # If we didn't find any, core_bases becomes an empty set.
    core_bases = set.union(*core_sets) if core_sets else set()

    # If the main CMV CSVs didn’t give us any sample IDs (e.g., missing files / empty outputs),
    # fall back to whatever we can detect from BLAST outputs and resistance tables.
    if not core_bases:
        if blast_top5 is not None and not blast_top5.empty and "base" in blast_top5.columns:
            core_bases |= set(blast_top5["base"].dropna().unique())
        if blast_hits_all is not None and not blast_hits_all.empty and "base" in blast_hits_all.columns:
            core_bases |= set(blast_hits_all["base"].dropna().unique())
        if resist_df is not None and not resist_df.empty and "base" in resist_df.columns:
            core_bases |= set(resist_df["base"].dropna().unique())

    # Sorted list of all discovered canonical sample IDs.
    all_bases = sorted(core_bases)

    def keep_bases(df):
        """
        Restrict a DataFrame to only rows belonging to the "core" set of sample IDs.

        Why this exists:
        - Different input sources may contain extra or unexpected sample IDs.
        - Filtering everything to core_bases keeps the app consistent so plots/tables
          don't include "orphan" samples that only appear in one file.

        Parameters
        ----------
        df : pd.DataFrame | None
            A DataFrame expected to contain a 'base' column.

        Returns
        -------
        pd.DataFrame | None
            - Original df if it is None/empty/doesn't have 'base'
            - Otherwise, a filtered copy containing only bases in core_bases
        """
        if df is None or df.empty or "base" not in df.columns:
            return df
        return df[df["base"].isin(core_bases)].copy()

    # Apply the "core base" restriction to every loaded table.
    # This prevents downstream sections from seeing mismatched sample lists.
    gene_df = keep_bases(gene_df)
    per_base_df = keep_bases(per_base_df)
    bam_df = keep_bases(bam_df)
    blast_top5 = keep_bases(blast_top5)
    blast_hits_all = keep_bases(blast_hits_all)
    align_df = keep_bases(align_df)
    vaf_df = keep_bases(vaf_df)
    resist_df = keep_bases(resist_df)

    # If we still have no samples, there is nothing meaningful to display.
    if not all_bases:
        st.warning("No samples found across inputs. Check your folder paths and file naming.")
        st.stop()

    # Optional debug panel: shows which sample IDs were detected from each source.
    # Useful when troubleshooting missing samples or naming mismatches.
    with st.sidebar.expander("Debug: bases detected by source", expanded=False):

        def names(df):
            """
            Return sorted unique base IDs from a DataFrame (or [] if not possible).

            This is purely for debugging and helps confirm that each input source
            contains the samples you expect.
            """
            return (
                sorted(df["base"].dropna().unique())
                if (df is not None and not df.empty and "base" in df.columns)
                else []
            )

        st.write(
            {
                "gene.csv": names(gene_df),
                "per_base.csv": names(per_base_df),
                "summary.csv": names(bam_df),
                "blast_top5.csv": names(blast_top5),
                "*.tsv (BLAST hits)": names(blast_hits_all),
                "BAM stats": names(align_df),
                "VAF tables": names(vaf_df),
                "Annotated VCFs": names(resist_df),
            }
        )

    # Sample selection widget:
    # - options = all detected bases
    # - default = first sample (so the app renders something immediately)
    with st.sidebar:
        bases_selected = st.multiselect(
            "Select sample(s)",
            options=all_bases,
            default=[all_bases[0]] if all_bases else [],
            help="Pick one or more canonical sample IDs.",
        )

    # If the user deselects everything, stop early to avoid errors in plots/tables.
    if not bases_selected:
        st.info("Select at least one base sample.")
        st.stop()

    def _keep_selected(df):
        """
        Filter a DataFrame down to the user-selected samples only.

        Notes
        -----
        - This assumes 'base' is the canonical sample column used throughout the app.
        - If df is missing or doesn't contain 'base', we return an empty DataFrame
          to keep downstream code simple (and avoid attribute errors).
        """
        return (
            df[df["base"].isin(bases_selected)].copy()
            if (df is not None and not df.empty and "base" in df.columns)
            else pd.DataFrame()
        )

    # Create filtered "working" tables used by the rest of the app for the chosen samples.
    f_align = _keep_selected(align_df)
    f_gene = _keep_selected(gene_df)

    # (CHANGE 2/3) only load per-base data when those sections are enabled
    #
    # per_base_df is often the largest file (many rows per genome position).
    # Loading it only when needed reduces startup time and makes the app feel snappier.
    if show_per_base_heatmaps or show_genome_wide:
        f_per = load_per_base_for_selected(csv_dir, tuple(bases_selected))
    else:
        f_per = pd.DataFrame()

    f_bam = _keep_selected(bam_df)
    f_vaf = _keep_selected(vaf_df)
    f_top5 = _keep_selected(blast_top5)
    f_res = _keep_selected(resist_df)

    # BLAST: restrict to CMV-ish titles and selected bases, then apply filters (NO E-VALUE)
    #
    # sel_hits_raw is the subset of BLAST hits that look like CMV/Merlin.
    # This is done using 'stitle' string matching (pragmatic filtering).
    sel_hits_raw = (
        blast_hits_all[
            blast_hits_all["stitle"].astype(str).str.contains(
                "Human herpesvirus 5|Human cytomegalovirus|Cytomegalovirus Merlin|NC_006273.2",
                case=False,
                na=False,
            )
        ].copy()
        if (blast_hits_all is not None and not blast_hits_all.empty)
        else pd.DataFrame()
    )

    # Restrict to only the user-selected bases (so plots/summaries don't include other samples).
    sel_hits_raw = (
        sel_hits_raw[sel_hits_raw["base"].isin(bases_selected)]
        if not sel_hits_raw.empty
        else sel_hits_raw
    )

    # Apply the BLAST filtering rules that your app uses (identity and alignment length).
    # This function is assumed to return BLAST hits that pass thresholds.
    sel_hits = filter_blast_merlin(
        sel_hits_raw,
        min_identity=float(min_identity),
        min_length=int(min_aln_len),
    )

    # ----------------------- 0) Analysis Summary (always shown) -----------------------
    #
    # This section creates a per-sample summary table:
    # - Alignment QC metrics (mapped %, properly paired %, duplicates %, singletons %)
    # - Coverage issues per gene
    # - Mean BLAST % identity for CMV-filtered hits
    # - Non-CMV species presence (from BLAST)
    # - CMV presence call using mapped_hq counts
    # - Mixed strain inference from mid-frequency variants
    # - Antiviral resistance reporting (3-state: detected / not detected / indeterminate)
    with st.expander("0. Analysis Summary", expanded=True):

        def _agg_align_by_base(f_align: pd.DataFrame) -> pd.DataFrame:
            """
            Aggregate alignment statistics per base (sample).

            We sum counts (e.g., mapped reads) and average insert size metrics.
            Then we derive percentages that are easier to interpret in QC.

            Returns
            -------
            pd.DataFrame
                One row per base, including:
                - total_reads, mapped_reads, duplicates, etc.
                - mapped_pct, properly_paired_pct, singletons_pct, duplicates_pct
            """
            if f_align is None or f_align.empty:
                return pd.DataFrame()

            # Group and aggregate numeric metrics by sample (base).
            grp = f_align.groupby("base", as_index=False).agg(
                {
                    "total_reads": "sum",
                    "mapped_reads": "sum",
                    "properly_paired_reads": "sum",
                    "paired_in_sequencing": "sum",
                    "singletons": "sum",
                    "duplicates": "sum",
                    "insert_size_mean": "mean",
                    "insert_size_sd": "mean",
                }
            )

            # Derived percentages. np.where avoids division by zero.
            grp["mapped_pct"] = np.where(
                grp["total_reads"] > 0,
                grp["mapped_reads"] / grp["total_reads"] * 100.0,
                np.nan,
            )
            grp["properly_paired_pct"] = np.where(
                grp["paired_in_sequencing"] > 0,
                grp["properly_paired_reads"] / grp["paired_in_sequencing"] * 100.0,
                np.nan,
            )
            grp["singletons_pct"] = np.where(
                grp["paired_in_sequencing"] > 0,
                grp["singletons"] / grp["paired_in_sequencing"] * 100.0,
                np.nan,
            )
            grp["duplicates_pct"] = np.where(
                grp["total_reads"] > 0,
                grp["duplicates"] / grp["total_reads"] * 100.0,
                np.nan,
            )
            return grp

        def _gene_depth_issues_by_base(
            f_gene: pd.DataFrame,
            th_mean_depth_min: float,
        ) -> dict:
            """
            Identify gene-level coverage issues per base.

            This looks at mean depth per gene and flags genes below a threshold.
            It returns a lightweight summary to display in the main table.

            Parameters
            ----------
            f_gene : pd.DataFrame
                Gene-level coverage table. Expected columns include 'base', 'gene', 'avg_depth'.
            th_mean_depth_min : float
                Minimum acceptable mean depth (×) per gene.

            Returns
            -------
            dict
                Mapping:
                    base -> (number_of_low_depth_genes, example_gene_names_list)
            """
            out = {}
            if f_gene is None or f_gene.empty:
                return out

            # Average depth per (base, gene). numeric_only=True avoids issues if avg_depth is read as string.
            g_agg = f_gene.groupby(["base", "gene"], as_index=False)["avg_depth"].mean(
                numeric_only=True
            )

            for b, gb in g_agg.groupby("base"):
                # Convert avg_depth safely, then filter below threshold.
                low = gb[pd.to_numeric(gb["avg_depth"], errors="coerce") < th_mean_depth_min]
                n = int(low.shape[0])
                names_ = list(low["gene"].astype(str).head(5))  # only show a few examples in summary
                out[b] = (n, names_)
            return out

        def _mean_identity_by_base(sel_hits_cmv: pd.DataFrame) -> pd.Series:
            """
            Compute mean BLAST % identity per base using CMV-filtered hits.

            Notes
            -----
            - Uses 'pident' as the % identity column.
            - Returns an empty Series if there are no hits.

            Returns
            -------
            pd.Series
                Index = base, value = mean % identity (float)
            """
            if (
                sel_hits_cmv is None
                or sel_hits_cmv.empty
                or "pident" not in sel_hits_cmv.columns
            ):
                return pd.Series(dtype=float)

            # Coerce to numeric, drop missing, then take mean.
            return sel_hits_cmv.groupby("base")["pident"].apply(
                lambda s: pd.to_numeric(s, errors="coerce").dropna().mean()
            )

        # Build summary output row-by-row to keep logic explicit and easy to follow.
        rows = []

        align_agg = _agg_align_by_base(f_align)
        depth_issues = _gene_depth_issues_by_base(f_gene, th_mean_depth_min)
        mean_ident = _mean_identity_by_base(sel_hits)

        # CMV mapped HQ reads per base (this is your CMV detection signal).
        # _agg_cmv_mapped_hq_by_base() is assumed to return columns:
        #   base, cmv_mapped_hq
        cmv_mapped_by_base = _agg_cmv_mapped_hq_by_base(f_bam)
        cmv_mapped_lookup = (
            {r["base"]: float(r["cmv_mapped_hq"]) for _, r in cmv_mapped_by_base.iterrows()}
            if not cmv_mapped_by_base.empty
            else {}
        )

        # Resistance details: build a lookup of base -> textual details.
        # This assumes 'details' holds your human-readable resistance call.
        res_map = {}
        if not f_res.empty:
            for _, r in f_res.iterrows():
                res_map[r["base"]] = r["details"]

        # Create one summary record per selected sample.
        for b in bases_selected:
            # Pull the alignment row for this base if available.
            arow = (
                align_agg[align_agg["base"] == b].iloc[0]
                if (not align_agg.empty and (align_agg["base"] == b).any())
                else None
            )

            # Extract key QC percentages, defaulting to NaN if unavailable.
            map_pct = float(arow["mapped_pct"]) if arow is not None and pd.notna(arow["mapped_pct"]) else np.nan
            prop_pct = float(arow["properly_paired_pct"]) if arow is not None and pd.notna(arow["properly_paired_pct"]) else np.nan
            dup_pct = float(arow["duplicates_pct"]) if arow is not None and pd.notna(arow["duplicates_pct"]) else np.nan
            sing_pct = float(arow["singletons_pct"]) if arow is not None and pd.notna(arow["singletons_pct"]) else np.nan

            # Alignment QC rules:
            # We build a list of "fail reasons" so we can show a clear summary.
            align_fail_reasons = []
            if not np.isnan(map_pct) and map_pct < th_map_pct_min:
                align_fail_reasons.append(f"mapped {map_pct:.1f}% < {th_map_pct_min}%")
            if not np.isnan(prop_pct) and prop_pct < th_proper_pct_min:
                align_fail_reasons.append(f"proper {prop_pct:.1f}% < {th_proper_pct_min}%")
            if not np.isnan(dup_pct) and dup_pct > th_dup_pct_max:
                align_fail_reasons.append(f"dup {dup_pct:.1f}% > {th_dup_pct_max}%")
            if not np.isnan(sing_pct) and sing_pct > th_singletons_max:
                align_fail_reasons.append(f"singletons {sing_pct:.1f}% > {th_singletons_max}%")

            align_status = "OK" if not align_fail_reasons else "Issues: " + " · ".join(align_fail_reasons)

            # Coverage QC based on gene-level depth.
            n_low, names_ = depth_issues.get(b, (0, []))
            cov_status = (
                "OK"
                if n_low == 0
                else f"{n_low} gene(s) < {th_mean_depth_min}×"
                + (f" (e.g., {', '.join(names_)})" if names_ else "")
            )

            # BLAST identity QC (mean % identity across CMV-filtered hits).
            mid = mean_ident.get(b, np.nan)
            id_status = (
                "OK"
                if (not np.isnan(mid) and mid >= th_identity_min)
                else ("No CMV hits" if np.isnan(mid) else f"Below threshold ({mid:.1f}% < {th_identity_min}%)")
            )
            id_text = f"{mid:.1f}%" if not np.isnan(mid) else "n/a"

            # Non-CMV detection using BLAST hits (separate helper function).
            # This is typically used to flag possible co-infections or contamination.
            non_cmv = (
                detect_non_cmv_species(blast_hits_all, b)
                if (blast_hits_all is not None and not blast_hits_all.empty)
                else False
            )
            non_cmv_status = "present" if non_cmv else "none"

            # CMV presence call using mapped HQ reads.
            # This is your primary detection output (Detected / Indeterminate / Not detected).
            cmv_mapped_hq = cmv_mapped_lookup.get(b, 0.0)
            cmv_status = cmv_presence_call_mapped_hq(
                cmv_mapped_hq,
                detected_threshold=float(th_detect_mapped),
                indeterminate_low=float(th_indeterminate_low),
                not_detected_high=float(th_not_detected_high),
            )
            presence_col = (
                f"{cmv_status} (mapped_hq={int(cmv_mapped_hq) if np.isfinite(cmv_mapped_hq) else 0})"
            )

            # Mixed strain inference based on VAF distribution (mid-frequency variants).
            # This is a heuristic: multiple strains often create clusters of mid-AF variants.
            vf_b = (
                f_vaf[f_vaf["base"] == b]
                if ("base" in f_vaf.columns and not f_vaf.empty)
                else pd.DataFrame()
            )
            mixed_verdict = infer_mixed_strains(
                vf_b,
                lo=th_mid_af_low,
                hi=th_mid_af_high,
                count_likely=th_mid_count_likely,
                frac_likely=th_mid_frac_likely,
                count_possible=th_mid_count_poss,
                frac_possible=th_mid_frac_poss,
            )

            # Antiviral resistance (3-state)
            # 1) If we have a resistance call from annotated VCFs, show it.
            # 2) Otherwise, decide whether we can *assess* resistance based on UL97/UL54 coverage.
            resist_text = res_map.get(b, "")
            if resist_text:
                resist_col = f"⚠️ {resist_text}"
            else:
                assessable, assess_detail = resistance_assessability_by_base(
                    f_gene, b, float(th_res_gene_depth_min)
                )
                genes_ok = int(assess_detail.get("genes_assessable", 0))
                ul97_md = assess_detail.get("UL97_mean_depth", np.nan)
                ul54_md = assess_detail.get("UL54_mean_depth", np.nan)

                # require at least N of {UL97, UL54} to meet depth threshold
                is_ok = genes_ok >= int(th_res_gene_min_count)

                if is_ok:
                    resist_col = "Not detected"
                else:
                    # Provide the mean depths to explain why it is indeterminate.
                    ul97_s = f"{ul97_md:.1f}×" if np.isfinite(ul97_md) else "n/a"
                    ul54_s = f"{ul54_md:.1f}×" if np.isfinite(ul54_md) else "n/a"
                    resist_col = (
                        "Indeterminate (insufficient UL97/UL54 coverage; "
                        f"UL97={ul97_s}, UL54={ul54_s})"
                    )

            # Add the per-base summary row.
            rows.append(
                {
                    "base": b,
                    "CMV detected": presence_col,
                    "Alignment": align_status if align_status == "OK" else f"⚠️ {align_status}",
                    "Coverage (mean per gene)": "OK" if cov_status == "OK" else f"⚠️ {cov_status}",
                    "CMV % identity (mean)": f"{id_status} ({id_text})" if id_status != "OK" else f"OK ({id_text})",
                    "Mixed strains": mixed_verdict,
                    "Non-CMV hits": non_cmv_status,
                    "Antiviral resistance": resist_col,
                }
            )

        # Convert rows into a DataFrame for display in the app.
        summary_df = pd.DataFrame(rows)
        st.dataframe(summary_df, use_container_width=True)

    # ----------------------- 1) Alignment Statistics (BWA-MEM) -----------------------
    # Shows the raw alignment stats table and a set of QC plots.
    if show_alignment:
        with st.expander("1. Alignment Statistics (BWA-MEM)", expanded=True):
            if not f_align.empty:
                st.dataframe(f_align.drop(columns=["base"]), use_container_width=True)

                # 4 columns of plots (if data are present).
                c1, c2, c3, c4 = st.columns(4)
                with c1:
                    fig = plot_alignment_stacked(f_align)
                    if fig:
                        st.plotly_chart(fig, use_container_width=True)
                with c2:
                    fig = plot_alignment_pct(f_align, "mapped_pct", "Alignment Rate (Mapped %)")
                    if fig:
                        st.plotly_chart(fig, use_container_width=True)
                with c3:
                    fig = plot_alignment_pct(f_align, "properly_paired_pct", "Properly Paired Reads (%)")
                    if fig:
                        st.plotly_chart(fig, use_container_width=True)
                with c4:
                    # Only plot insert size if we actually have values.
                    if "insert_size_mean" in f_align.columns and f_align["insert_size_mean"].notna().any():
                        fig = px.bar(
                            f_align,
                            x="sample",
                            y="insert_size_mean",
                            error_y="insert_size_sd",
                            title="Insert Size (mean ± SD)",
                        )
                        st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No alignment stats found for selected samples.")

    # ----------------------- 2. Average Depth per Gene -----------------------
    # Displays gene-level coverage (mean depth per gene) and a plot summarising it.
    if show_depth_per_gene:
        with st.expander("2. Average Depth per Gene", expanded=True):
            if not f_gene.empty:
                st.dataframe(f_gene.drop(columns=["base"]), use_container_width=True)
                fig = plot_avg_depth(f_gene)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No gene-level coverage data for selected samples.")

    # ----------------------- 3. BLAST Summary, % Identity -----------------------
    # BLAST-derived summaries/plots (note: this implementation intentionally avoids E-values).
    if show_blast:
        with st.expander("3. BLAST Summary & % Identity (no E-values)", expanded=True):
            if show_blast_bar and not f_top5.empty:
                st.subheader("Top BLAST species per sample (selected)")
                st.dataframe(f_top5.drop(columns=["base"]), use_container_width=True)

                # Bar plot of top species hit counts by sample.
                if "species" in f_top5.columns and "total_hits" in f_top5.columns:
                    fig = px.bar(
                        f_top5,
                        x="sample",
                        y="total_hits",
                        color="species",
                        title="Top BLAST Species per Sample and Read",
                    )
                    fig.update_xaxes(tickangle=45)
                    st.plotly_chart(fig, use_container_width=True)
            elif show_blast_bar and f_top5.empty:
                st.info("Top BLAST species table not available.")

            # Genome-wide BLAST identity heatmap.
            if show_blast_heatmap and not sel_hits.empty:
                st.subheader("BLAST % identity across genome (selected, filtered)")
                fig = plot_blast_full_heatmap(sel_hits, bin_size=100)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            elif show_blast_heatmap and (sel_hits is None or sel_hits.empty):
                st.info("No CMV-filtered BLAST hits for genome heatmap.")

            # BLAST identity by predefined conserved regions.
            if show_blast_regions and not sel_hits.empty:
                st.subheader("BLAST % identity by conserved region (selected, filtered)")
                tabs = st.tabs([r["name"] for r in CONSERVE_REGIONS])
                for tab, region in zip(tabs, CONSERVE_REGIONS):
                    with tab:
                        fig = plot_blast_region_heatmap(sel_hits, region=region)
                        if fig:
                            st.plotly_chart(fig, use_container_width=True)
            elif show_blast_regions and (sel_hits is None or sel_hits.empty):
                st.info("No CMV-filtered BLAST hits for region heatmaps.")

    # ----------------------- 4. Per-base Coverage Heatmaps -----------------------
    # Per-base depth heatmaps across conserved regions (requires f_per).
    if show_per_base_heatmaps:
        with st.expander("4. Per-base Coverage Heatmaps", expanded=False):
            if not f_per.empty:
                tabs = st.tabs([r["name"] for r in CONSERVE_REGIONS])
                for tab, region in zip(tabs, CONSERVE_REGIONS):
                    with tab:
                        fig = plot_per_base_heatmap(f_per, gene=region["name"], bin_size=None)
                        if fig:
                            st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No per-base depth data for selected samples.")

    # ----------------------- 5. Genome-wide Coverage -----------------------
    # Whole-genome depth plot (requires f_per).
    if show_genome_wide:
        with st.expander("5. Genome-wide Coverage", expanded=False):
            if not f_per.empty:
                fig = plot_genome_coverage(f_per)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No per-base depth data for genome-wide view (selected).")

    # ----------------------- 6. FastQC Reports -----------------------
    # Lists FastQC outputs and allows previewing HTML reports directly in Streamlit.
    if show_fastqc:
        with st.expander("6. FastQC Reports", expanded=False):
            if fastqc_dir and os.path.isdir(fastqc_dir):
                selected_bases_set = set(bases_selected)

                def _is_selected_fastqc(path: str) -> bool:
                    """
                    Decide whether a FastQC file belongs to a selected base.

                    FastQC outputs are typically named like:
                        <sample>_fastqc.html
                        <sample>_fastqc.zip

                    We strip the suffix, then convert to canonical sample naming using
                    canonical_sample(), then check it is in the selected set.
                    """
                    fname = os.path.basename(path)
                    sample = fname.replace("_fastqc.html", "").replace("_fastqc.zip", "")
                    return canonical_sample(sample) in selected_bases_set

                html_files_all = sorted(glob.glob(os.path.join(fastqc_dir, "*_fastqc.html")))
                zip_files_all = sorted(glob.glob(os.path.join(fastqc_dir, "*_fastqc.zip")))
                html_files = [p for p in html_files_all if _is_selected_fastqc(p)]
                zip_files = [p for p in zip_files_all if _is_selected_fastqc(p)]

                # Build a module-level summary across selected FastQC ZIPs.
                all_sum = []
                for z in zip_files:
                    try:
                        sample = os.path.basename(z).replace("_fastqc.zip", "")
                        with zipfile.ZipFile(z) as zf:
                            inner = [n for n in zf.namelist() if n.endswith("/summary.txt")]
                            if inner:
                                with zf.open(inner[0]) as fh:
                                    for line in io.TextIOWrapper(
                                        fh, encoding="utf-8", errors="ignore"
                                    ):
                                        parts = line.strip().split("\t")
                                        if len(parts) >= 3:
                                            status, module, _filename = parts[0], parts[1], parts[2]
                                            all_sum.append(
                                                {"sample": sample, "module": module, "status": status}
                                            )
                    except Exception as e:
                        st.warning(f"Failed to read {z}: {e}")

                # Display a pivoted table (sample x module) of PASS/WARN/FAIL summaries.
                if all_sum:
                    df_sum = pd.DataFrame(all_sum)
                    df_sum["base"] = df_sum["sample"].apply(canonical_sample)
                    df_sum = df_sum[df_sum["base"].isin(bases_selected)].drop(columns=["base"])
                    st.subheader("FastQC module summary (selected)")
                    pivot = df_sum.pivot_table(
                        index=["sample"],
                        columns="module",
                        values="status",
                        aggfunc=lambda x: ", ".join(sorted(set(x))),
                    )
                    st.dataframe(pivot.fillna("N/A"), use_container_width=True)

                # Allow user to preview either HTML files directly or HTML inside ZIP files.
                reports = html_files + zip_files
                if reports:
                    choice = st.selectbox("Preview a FastQC report", options=["(select)"] + reports)
                    if choice and choice != "(select)":
                        st.markdown(f"**Preview:** {os.path.basename(choice)}")
                        try:
                            if choice.endswith(".zip"):
                                with zipfile.ZipFile(choice) as zf:
                                    html_inside = [f for f in zf.namelist() if f.endswith("fastqc_report.html")]
                                    if html_inside:
                                        with zf.open(html_inside[0]) as f:
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
    # Displays VAF table and two plots:
    # - AF across genome (scatter)
    # - AF distribution (histogram)
    if show_variants:
        with st.expander("7. Variants (LoFreq VAF)", expanded=True):
            if f_vaf is None or f_vaf.empty:
                st.info("No VAF tables loaded for selected samples.")
            else:
                st.subheader("VAF table (selected samples)")
                st.dataframe(
                    f_vaf[
                        [c for c in ["sample", "CHROM", "POS", "REF", "ALT", "DP", "AD_ALT", "AF", "base"] if c in f_vaf.columns]
                    ]
                    .sort_values(["sample", "POS"])
                    .reset_index(drop=True)
                    .drop(columns=["base"]),
                    use_container_width=True,
                )

                c1, c2 = st.columns(2)
                with c1:
                    # Scatter plot of AF by genome position, per sample.
                    fig = go.Figure()
                    for sample, sdf in f_vaf.groupby("sample"):
                        sdf = sdf.sort_values("POS")
                        fig.add_trace(
                            go.Scatter(
                                x=sdf["POS"],
                                y=sdf["AF"],
                                mode="markers",
                                name=sample,
                                hovertemplate="POS=%{x}<br>AF=%{y:.3f}<extra>",
                            )
                        )
                    fig.update_layout(
                        title="Variant Allele Frequency (AF) across genome",
                        xaxis_title="Position (bp)",
                        yaxis_title="AF (0-1)",
                        yaxis=dict(range=[0, 1]),
                    )
                    st.plotly_chart(fig, use_container_width=True)

                with c2:
                    # Histogram to show how AF values are distributed per sample.
                    fig = px.histogram(f_vaf, x="AF", color="sample", nbins=40, title="Distribution of AF values")
                    fig.update_xaxes(range=[0, 1])
                    st.plotly_chart(fig, use_container_width=True)

# If the user has entered a CSV directory but hasn't clicked "Load / Refresh" yet,
# we show instructions rather than attempting to load anything.
elif csv_dir and not st.session_state.get("_load_committed", False):
    st.info("Set your directories, then click **Load / Refresh** in the sidebar to load data.")
else:
    # No CSV directory provided, so we can't proceed.
    st.info("Enter at least a CSV directory to load data.")
