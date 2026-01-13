#!/usr/bin/env python3
"""
CMV pipeline GUI runner/

Functionality overview: 
Desktop GUI that runs CMV sequencing pipeline end-to-end.
---------------------------
- Pick input folders/files (FASTQs, references, BLAST DB, etc.) using Browse buttons.
- Run selected pipeline steps.
- Stream logs live into the GUI window.
- Optionally compress FASTQ outputs.
- Generate the alignment-metrics CSV that the Streamlit report expects (summary.csv).
- Launch the Streamlit interactive report on a chosen port
- Open FastQC summary.html and launch IGV for quick review
"""

import os
import re
import signal
import subprocess
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
import webbrowser
import glob
import matplotlib

matplotlib.use('Agg')  # Use non-GUI backend for matplotlib
import shlex
import datetime
import shutil
import tempfile
import tarfile
import zipfile
import gzip
import csv
from pathlib import Path

from extract_cmv_bam_metrics import run_bam_metrics
from summarize_blast_hits import summarize_blast_hits

# ---------- archive/fastq constants ----------
# These are used to decide how a user-selected path should be handled:
# - if it's an archive -> unpack it
# - if it's a gzipped FASTA -> gunzip it
# - if it's a FASTQ -> accept as input
ARCHIVE_EXTS = ('.zip', '.tar', '.tar.gz', '.tgz', '.tar.bz2', '.tbz2', '.tar.xz', '.txz')
GZIP_FASTA_EXTS = ('.fa.gz', '.fasta.gz', '.fna.gz')
FASTQ_EXTS = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

# -------------------- sample name standardisation --------------------
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
    r"(?:_indelqual)$",
    r"(?:_realigned)$",
    r"(?:\.bam)$",
]


def canonical_sample(name: str) -> str:
    """
    Convert a filename-ish sample string into a stable "base" sample name.

    This repeatedly strips common suffix patterns used by sequencing/pipeline outputs
    (R1/R2 tags, sorted/dedup suffixes, ".bam", etc.) until nothing else can be removed.

    Parameters
    ----------
    name : str
        Input name (often a BAM filename without the directory).

    Returns
    -------
    str
        Standard sample name (used as the 'base' key in summary outputs).
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


# -------------------- samtools-based CMV alignment metrics (HQ mapped + breadth) --------------------
def _run_cmd(cmd: list[str]) -> str:
    """
    Run a command and return stdout as text.

     subprocess.run is wrapped so we can:
    - capture stdout/stderr
    - raise a clear error if the command fails (return code != 0)

    Parameters
    ----------
    cmd : list[str]
        Command and arguments, e.g. ["samtools", "view", "-c", "file.bam"].

    Returns
    -------
    str
        Standard output (stdout) from the command.

    Raises
    ------
    RuntimeError
        If the command returns a non-zero exit code.
    """
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n  {' '.join(cmd)}\n\nSTDERR:\n{p.stderr}")
    return p.stdout


def samtools_count(bam: str, mapped_only: bool = False, mapq: int | None = None) -> int:
    """
    Count reads in a BAM using `samtools view -c`.

    Parameters
    ----------
    bam : str
        Path to BAM.
    mapped_only : bool, default=False
        If True, excludes unmapped reads (samtools flag -F 4).
    mapq : int | None, default=None
        If provided, only count reads with mapping quality >= mapq (samtools -q).

    Returns
    -------
    int
        Read count matching the filters.
    """
    cmd = ["samtools", "view", "-c"]
    if mapped_only:
        cmd += ["-F", "4"]  # exclude unmapped reads
    if mapq is not None:
        cmd += ["-q", str(int(mapq))]  # MAPQ filter
    cmd += [bam]
    out = _run_cmd(cmd).strip()
    return int(out) if out else 0


def samtools_idxstats_lengths(bam: str) -> int:
    """
    Sum reference contig lengths using `samtools idxstats`.

    This is used to estimate the reference length for breadth calculations.

    Parameters
    ----------
    bam : str
        Path to BAM (must be indexed for idxstats to work reliably).

    Returns
    -------
    int
        Total reference length across contigs (excluding the '*' line).
    """
    out = _run_cmd(["samtools", "idxstats", bam])
    genome_len = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        contig, length = parts[0], parts[1]
        if contig == "*":
            continue
        try:
            genome_len += int(length)
        except ValueError:
            continue
    return genome_len


def breadth_of_coverage(
    bam: str,
    mapq: int,
    depth_min: int = 1,
    method: str = "coverage",
) -> tuple[float, int, int]:
    """
    Calculate breadth of coverage for a BAM.

    Parameters
    ----------
    bam : str
        Path to BAM.
    mapq : int
        Minimum mapping quality for included reads.
    depth_min : int, default=1
        Minimum depth to call a position "covered".
    method : str, default="coverage"
        Preferred method:
        - "coverage": uses `samtools coverage` (usually faster)
        - falls back to "depth": uses `samtools depth -a` (slower, but robust)

    Returns
    -------
    tuple[float, int, int]
        (breadth, covered_bases, ref_bases)
        - breadth may be NaN if reference length cannot be determined
    """
    # First get the reference size from idxstats
    ref_bases = samtools_idxstats_lengths(bam)
    if ref_bases <= 0:
        return float("nan"), 0, 0

    # Preferred: samtools coverage
    if method == "coverage":
        try:
            out = _run_cmd(["samtools", "coverage", "-q", str(int(mapq)), bam])
            covbases_sum = 0
            len_sum = 0
            for line in out.splitlines():
                if not line or line.startswith("#") or line.lower().startswith("rname"):
                    continue
                parts = line.split()
                # Expected columns include: rname start end numreads covbases coverage meandepth ...
                if len(parts) < 6:
                    continue
                rname = parts[0]
                if rname == "*":
                    continue
                start = int(parts[1])
                end = int(parts[2])
                covbases = int(parts[4])
                length = end - start + 1
                covbases_sum += covbases
                len_sum += length
            if len_sum > 0:
                breadth = covbases_sum / len_sum
                return float(breadth), int(covbases_sum), int(len_sum)
        except Exception:
            # If samtools coverage fails (older version, indexing issues, etc.),
            # fall back to samtools depth which is slower but more universal.
            method = "depth"

    # Fallback: samtools depth -a (slower)
    out = _run_cmd(["samtools", "depth", "-a", "-q", str(int(mapq)), bam])
    covered = 0
    total = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            d = int(parts[2])
        except ValueError:
            continue
        total += 1
        if d >= depth_min:
            covered += 1

    if total == 0:
        return float("nan"), 0, 0
    return float(covered / total), int(covered), int(total)


def write_cmv_alignment_metrics_csv(
    bam_dir: str,
    out_dir: str,
    mapq: int = 30,
    depth_min: int = 1,
    breadth_method: str = "coverage",
    log_func=None,
) -> bool:
    """
    Write the alignment-based summary CSV files the Streamlit report expects.

    This function scans a folder of CMV-aligned BAMs and produces a per-sample summary with:
    - total reads in BAM (all)
    - HQ mapped reads (MAPQ >= mapq, excluding unmapped)
    - RPM using HQ mapped reads
    - breadth of coverage (with MAPQ filter and depth threshold)

    Output files (same content):
      - <out_dir>/summary.csv
      - <out_dir>/cmv_alignment_summary.csv

    Parameters
    ----------
    bam_dir : str
        Directory containing BAM files.
    out_dir : str
        Output directory to write summary files into.
    mapq : int, default=30
        Mapping quality threshold used for HQ mapped reads and breadth calculation.
    depth_min : int, default=1
        Minimum depth to consider a position covered.
    breadth_method : str, default="coverage"
        "coverage" preferred; "depth" fallback.
    log_func : callable | None
        Optional logging callback (e.g. GUI logger). If None, logging is skipped.

    Returns
    -------
    bool
        True on success.

    Raises
    ------
    RuntimeError
        If samtools is missing or no BAMs are found or any BAM processing fails.
    """
    def log(msg: str):
        # Local wrapper so this function can run both in GUI and non-GUI contexts.
        if log_func:
            log_func(msg)

    # samtools is required for view/idxstats/coverage/depth.
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found in PATH (required to compute HQ-mapped metrics).")

    os.makedirs(out_dir, exist_ok=True)

    # Collect BAMs (ignore BAM index files)
    bams = sorted([p for p in glob.glob(os.path.join(bam_dir, "*.bam")) if not p.endswith(".bai")])
    if not bams:
        raise RuntimeError(f"No .bam files found in: {bam_dir}")

    rows = []
    for bam in bams:
        # Derive the sample key ('base') from the BAM name
        sample = os.path.basename(bam)
        sample_noext = re.sub(r"\.bam$", "", sample, flags=re.IGNORECASE)
        base = canonical_sample(sample_noext)

        try:
            # total_reads_bam: all reads in BAM (mapped + unmapped)
            total_reads_bam = samtools_count(bam, mapped_only=False, mapq=None)

            # mapped_hq: mapped reads with MAPQ >= threshold
            mapped_hq = samtools_count(bam, mapped_only=True, mapq=mapq)

            # RPM: reads-per-million (normalises by total reads in this BAM)
            rpm_hq = (mapped_hq / total_reads_bam * 1e6) if total_reads_bam else 0.0

            # Breadth of coverage: fraction of positions covered
            breadth, cov_bases, ref_bases = breadth_of_coverage(
                bam, mapq=mapq, depth_min=depth_min, method=breadth_method
            )

            rows.append({
                "base": base,
                "bam": bam,
                "total_reads_in_cmv_bam": int(total_reads_bam),
                "cmv_mapped_reads_mapq": int(mapped_hq),
                "cmv_rpm_mapq": float(rpm_hq),
                "cmv_breadth_mapq": float(breadth),
                "cmv_cov_bases_mapq": int(cov_bases),
                "cmv_ref_bases": int(ref_bases),
            })
            log(f"[metrics] {base}: total={total_reads_bam} mapped_hq(MAPQ≥{mapq})={mapped_hq} breadth={breadth:.4g}")
        except Exception as e:
            # Fail early because the summary file is only useful if all samples are correct.
            log(f"[metrics] ❌ Failed for {bam}: {e}")
            raise

    # Stable column order (what your app expects)
    fieldnames = [
        "base",
        "bam",
        "total_reads_in_cmv_bam",
        "cmv_mapped_reads_mapq",
        "cmv_rpm_mapq",
        "cmv_breadth_mapq",
        "cmv_cov_bases_mapq",
        "cmv_ref_bases",
    ]

    summary_path = os.path.join(out_dir, "summary.csv")
    safety_path = os.path.join(out_dir, "cmv_alignment_summary.csv")

    # Write both files so you always have a safe copy if something overwrites summary.csv later.
    for path in (summary_path, safety_path):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            # Sort by base so the output is stable/repeatable across runs
            for r in sorted(rows, key=lambda x: x["base"]):
                w.writerow(r)

    log(f"[metrics] ✅ Wrote alignment metrics CSV: {summary_path}")
    log(f"[metrics] ✅ Wrote alignment metrics CSV (copy): {safety_path}")
    return True


# ------------------------------- UI blocks -------------------------------

class CollapsibleStep(tk.Frame):
    """Collapsible block for pipeline."""
    def __init__(self, parent, step_name, labels, browse_types=None, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.step_name = step_name
        self.labels = labels
        self.browse_types = browse_types or {}

        # Each step is enabled/disabled via a checkbox (BooleanVar drives the UI state).
        self.var = tk.BooleanVar()

        # Step header checkbox
        self.checkbox = tk.Checkbutton(
            self,
            text=step_name.replace('_', ' ').title(),
            variable=self.var,
            command=self.toggle,
            font=("Arial", 11, "bold"),
        )
        self.checkbox.pack(anchor='w', pady=2)

        # Contents are placed in a child frame which we pack/unpack to collapse
        self.content_frame = tk.Frame(self)
        self.entries = {}

        # Build one label + entry row per parameter for this step
        for label in labels:
            row = tk.Frame(self.content_frame)
            lbl = tk.Label(row, text=f"{label}:", width=22, anchor='e')
            ent = tk.Entry(row, width=54)

            # Decide if this label should Browse for a directory/file or be plain text
            browse_type = self.browse_types.get(label.lower(), 'dir')
            if browse_type != 'text':
                btn = tk.Button(row, text="Browse", command=lambda l=label: self.browse_path(l))
                btn.pack(side='right', padx=5)

            lbl.pack(side='left')
            ent.pack(side='left', fill='x', expand=True, padx=5)
            row.pack(fill='x', pady=2)
            self.entries[label] = ent

        # Start collapsed by default
        self.content_frame.pack_forget()

    def toggle(self):
        """Show/hide the step input fields depending on the checkbox state."""
        if self.var.get():
            self.content_frame.pack(fill='x', padx=25, pady=5)
        else:
            self.content_frame.pack_forget()

    def browse_path(self, label):
        """
        Open file/directory dialogs based on the field type.

        Behaviour:
        - 'dir': askdirectory first; if cancelled, offer archive selection
        - 'file': askopenfilename with file-type filters based on label keywords
        - 'text': no browse (handled by not creating a button)
        """
        bt = self.browse_types.get(label.lower(), 'dir')
        lbl = label.lower()
        path = ""

        if bt == 'dir':
            # Directory chooser; if user cancels, allow selecting an archive instead
            path = filedialog.askdirectory(title=f"Select {label}")
            if not path:
                path = filedialog.askopenfilename(
                    title=f"Or select an archived {label}",
                    filetypes=[("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                               ("All files", "*.*")]
                )

        elif bt == 'file':
            # File chooser with filters tuned to common bioinformatics file types
            if "fastq" in lbl:
                filetypes = [("FASTQ files", "*.fastq *.fq *.fastq.gz *.fq.gz"),
                             ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            elif "fasta" in lbl or "reference" in lbl or "protein" in lbl:
                filetypes = [("FASTA files", "*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz"),
                             ("All files", "*.*")]
            elif "gff3" in lbl or "gff" in lbl:
                filetypes = [("GFF3 files", "*.gff3 *.gff"), ("All files", "*.*")]
            elif "bam" in lbl:
                filetypes = [("BAM files", "*.bam"),
                             ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            elif "tsv" in lbl:
                filetypes = [("TSV files", "*.tsv *.tsv.gz"),
                             ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            else:
                filetypes = [("All files", "*.*")]

            path = filedialog.askopenfilename(title=f"Select {label}", filetypes=filetypes)

        else:
            # Unknown browse type; do nothing
            return

        # Write chosen path into the entry box
        if path:
            self.entries[label].delete(0, tk.END)
            self.entries[label].insert(0, path)

    def get_values(self):
        """Return current entry values (trimmed) as a dict keyed by the labels."""
        return {label: ent.get().strip() for label, ent in self.entries.items()}


# ------------------------------- Main App -------------------------------

class PipelineGUI(tk.Tk):
    """
    Main Tkinter GUI application for running the CMV pipeline.

    The pipeline uses checkboxes (include/exclude each step) and a fixed step order.
    Execution happens in a background thread so the GUI stays responsive.

    (Per your note: unchanged; only snpEff bits added as optional)
    """
    def __init__(self):
        super().__init__()
        self.title("CMV Bioinformatics Pipeline GUI")
        self.geometry("980x980")

        # Holds the Streamlit process to avoid launching duplicates and terminate on close
        self._streamlit_proc = None

        # Output handling options
        self.compress_fastq = tk.BooleanVar(value=False)
        self.delete_uncompressed_fastq = tk.BooleanVar(value=False)
        self.use_pigz = tk.BooleanVar(value=True)
        self.tar_blast_output = tk.BooleanVar(value=False)

        # Track temporary staging directories created during archive extraction / gunzip steps
        self._temps = []

        # Single separate "step" for the Streamlit report launcher UI block
        self.steps_info = {"interactive_report_streamlit": ["Path to app.py", "Port (e.g., 8501)"]}

        # How each label should be browsed (file/dir/text)
        self.browse_types = {
            "path to app.py": "file",
            "port (e.g., 8501)": "text",
            "root output directory": "dir",
            "initial fastq directory": "dir",
            "human reference genome": "dir",
            "cmv reference genome": "dir",
            "blast db directory": "dir",
            "quality threshold": "text",
            "threads/cores": "text",
            # NEW optional snpEff inputs
            "gff3 annotation (for snpeff)": "file",
            "protein fasta (optional, snpeff)": "file",
        }

        self.steps = {}          # holds step widgets (like the Streamlit launcher block)
        self.one_click_vars = {} # holds tk.BooleanVar per pipeline step checkbox

        # Execution order matters; this list defines the pipeline sequence in the GUI
        self.one_click_step_order = [
            "split_interleaved",
            "filtering",
            "index_human",
            "align_human",
            "unaligned_reads",
            "index_cmv",
            "align_cmv",
            "extract_cmv_metrics",
            "blast_summary",
            "summarize_blast_hits",
            "variant_calling_vaf",
            "annotate_vcf_snpeff", 
        ]

        self.build_ui()
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    # ------------------------------ helpers: temps & archives ------------------------------
    def _register_temp(self, p: str) -> str:
        """Keep track of temporary paths to clean them up on exit."""
        self._temps.append(p)
        return p

    def _cleanup_temps(self):
        """Remove any temporary staging files/folders created during this GUI run."""
        for p in getattr(self, "_temps", []):
            try:
                if os.path.isdir(p):
                    shutil.rmtree(p, ignore_errors=True)
                elif os.path.isfile(p):
                    os.remove(p)
            except Exception:
                # Best-effort cleanup; do not crash on close
                pass
        self._temps = []

    def _is_archive(self, p: str) -> bool:
        """Return True if a path points to a file with a supported archive extension."""
        return bool(p) and os.path.isfile(p) and p.lower().endswith(ARCHIVE_EXTS)

    def _unpack_archive(self, archive_path: str, name_hint: str) -> str:
        """
        Unpack an archive into a temp directory and return a usable directory path.

        If the archive contains a single top-level folder, return that folder.
        Otherwise return the extraction directory itself.
        """
        dest = self._register_temp(tempfile.mkdtemp(prefix=f"cmv_{name_hint}_"))
        try:
            shutil.unpack_archive(archive_path, dest)
        except Exception as e:
            raise RuntimeError(f"Failed to unpack {archive_path}: {e}")
        entries = [os.path.join(dest, x) for x in os.listdir(dest)]
        dirs = [d for d in entries if os.path.isdir(d)]
        return dirs[0] if len(dirs) == 1 else dest

    def _gunzip_file(self, gz_path: str, name_hint: str) -> str:
        """Gunzip a single .gz file to a temp folder and return the decompressed file path."""
        out_dir = self._register_temp(tempfile.mkdtemp(prefix=f"cmv_{name_hint}_"))
        out_path = os.path.join(out_dir, os.path.basename(gz_path)[:-3])
        with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return out_path

    def _stage_single_file(self, src_file: str, name_hint: str) -> str:
        """
        Stage a single file into a temp folder.

        """
        d = self._register_temp(tempfile.mkdtemp(prefix=f"cmv_{name_hint}_"))
        dst = os.path.join(d, os.path.basename(src_file))
        try:
            os.symlink(os.path.abspath(src_file), dst)
        except Exception:
            shutil.copy2(src_file, dst)
        return d

    # ------------------------------ helpers: discovery ------------------------------
    def _ensure_dir(self, p: str):
        """Create a directory if needed and return the path."""
        os.makedirs(p, exist_ok=True)
        return p

    def _has_fastq(self, path: str) -> bool:
        """Check whether a directory contains FASTQ files (gz or not)."""
        pats = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
        return any(glob.glob(os.path.join(path, p)) for p in pats)

    def _auto_fastq_dir(self, base: str) -> str:
        """
        Try to locate FASTQs beneath a starting folder.

        """
        if self._has_fastq(base):
            return base

        for name in ("filtered_reads", "split", "fastq", "fastqs", "reads", "unaligned_reads", "split_interleaved"):
            cand = os.path.join(base, name)
            if os.path.isdir(cand) and self._has_fastq(cand):
                return cand

        try:
            for e in os.scandir(base):
                if e.is_dir() and self._has_fastq(e.path):
                    return e.path
        except Exception:
            pass

        return base

    def _assert_fastqs(self, d: str, stepname: str):
        """Raise an error early if a step expects FASTQs but none are found."""
        if not self._has_fastq(d):
            raise RuntimeError(f"{stepname}: no FASTQ files found in {d}.")

    def _ref_dir_and_file(self, ref_path: str):
        """
        Allow alignment scripts to accept either:
        - a folder containing reference + index
        - or a direct FASTA path
        """
        if os.path.isdir(ref_path):
            return ref_path, None
        return os.path.dirname(ref_path), ref_path

    def _find_bam_dir_from_alignment(self, align_folder: str) -> str:
        """
        Locate BAMs produced by alignment steps.

        Primary expectation: <align_folder>/BAM/*.bam
        Fallback: search recursively.
        """
        candidate = os.path.join(align_folder, "BAM")
        if os.path.isdir(candidate) and glob.glob(os.path.join(candidate, "*.bam")):
            return candidate
        hits = glob.glob(os.path.join(align_folder, "**", "*.bam"), recursive=True)
        if hits:
            return os.path.dirname(hits[0])
        raise RuntimeError(f"No BAM files found under {align_folder}")

    def _pick_fastq_from_unaligned_out(self, out_dir: str) -> str:
        """
        Unaligned reads step may output FASTQs in different layouts.
        This tries common locations to find them.
        """
        if self._has_fastq(out_dir):
            return out_dir
        nested = os.path.join(out_dir, "unaligned_reads")
        if os.path.isdir(nested) and self._has_fastq(nested):
            return nested
        return self._auto_fastq_dir(out_dir)

    def _symlink_or_copy_fastqs(self, src_dir: str, dst_dir: str):
        """
        Stage FASTQs into a working directory.

        This avoids modifying the original input folder and standardises the layout.
        """
        self._ensure_dir(dst_dir)
        patterns = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
        files = []
        for pat in patterns:
            files.extend(glob.glob(os.path.join(src_dir, pat)))
        if not files:
            raise RuntimeError(f"split_interleaved: no FASTQs to stage from {src_dir}")

        for f in files:
            dst = os.path.join(dst_dir, os.path.basename(f))
            if os.path.exists(dst):
                continue
            try:
                os.symlink(os.path.abspath(f), dst)
            except Exception:
                shutil.copy2(f, dst)

    # ------------------------------ helpers: input normalizers ------------------------------
    def _normalize_fastq_input(self, path: str) -> str:
        """
        Normalise FASTQ input path into a directory suitable for pipeline scripts.

        Supported inputs:
        - directory of FASTQs
        - archive containing FASTQs (unpacked to temp)
        - single FASTQ file (staged to temp)
        """
        if os.path.isdir(path):
            return self._auto_fastq_dir(path)
        if self._is_archive(path):
            return self._auto_fastq_dir(self._unpack_archive(path, "fastq"))
        if any(path.lower().endswith(ext) for ext in FASTQ_EXTS):
            return self._stage_single_file(path, "fastq_single")
        raise RuntimeError(f"FASTQ input not found/unsupported: {path}")

    def _normalize_bam_dir(self, path: str) -> str:
        """
        Normalise BAM input into a directory that contains BAMs.

        Supported inputs:
        - directory containing BAMs
        - archive containing BAMs (unpacked to temp)
        - single .bam file (staged to temp)
        """
        if os.path.isdir(path):
            if glob.glob(os.path.join(path, "*.bam")) or glob.glob(os.path.join(path, "**", "*.bam"), recursive=True):
                return path
            raise RuntimeError(f"No BAM files found in directory {path}")
        if self._is_archive(path):
            extracted = self._unpack_archive(path, "bam")
            hits = glob.glob(os.path.join(extracted, "**", "*.bam"), recursive=True)
            if not hits:
                raise RuntimeError(f"No BAM files found inside archive {path}")
            return os.path.dirname(hits[0])
        if path.lower().endswith(".bam"):
            return self._stage_single_file(path, "bam_single")
        raise RuntimeError(f"BAM input not found/unsupported: {path}")

    def _normalize_blast_db(self, path: str) -> str:
        """BLAST DB must be a folder or an archive that unpacks to a folder."""
        if os.path.isdir(path):
            return path
        if self._is_archive(path):
            return self._unpack_archive(path, "blastdb")
        raise RuntimeError(f"BLAST DB path must be a folder or an archive: {path}")

    def _normalize_tsv_dir(self, path: str) -> str:
        """
        Normalise BLAST TSV input into a directory.

        Supported:
        - folder of TSVs
        - archive of TSVs
        - single .tsv or .tsv.gz (staged to temp)
        """
        if os.path.isdir(path):
            return path
        if self._is_archive(path):
            return self._unpack_archive(path, "tsv")
        if path.lower().endswith((".tsv", ".tsv.gz")):
            return self._stage_single_file(path, "tsv_single")
        raise RuntimeError(f"BLAST TSV input not found/unsupported: {path}")

    def _normalize_reference_for_alignment(self, ref: str) -> str:
        """
        Normalise a reference path for alignment steps.

        Supported:
        - reference folder (contains index etc.)
        - archive -> unpack to temp
        - plain FASTA file
        - gzipped FASTA -> gunzip to temp and use the decompressed FASTA
        """
        if os.path.isdir(ref):
            return ref
        if self._is_archive(ref):
            return self._unpack_archive(ref, "ref")
        low = ref.lower()
        if low.endswith((".fa", ".fasta", ".fna")):
            return ref
        if low.endswith(GZIP_FASTA_EXTS):
            return self._gunzip_file(ref, "ref_gunzip")
        raise RuntimeError(f"Reference not found/unsupported: {ref}")

    def _ref_fasta_file(self, ref_input: str) -> str:
        """
        Resolve a reference input into an actual FASTA file path.

        This is used for VAF calling and snpEff DB building, which typically want a FASTA file.
        """
        ref_path = ref_input
        if self._is_archive(ref_input):
            ref_path = self._unpack_archive(ref_input, "ref_vaf")

        if os.path.isdir(ref_path):
            for pat in ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"):
                hits = glob.glob(os.path.join(ref_path, pat))
                if hits:
                    ref_path = hits[0]
                    break

        low = ref_path.lower()
        if low.endswith((".fa", ".fasta", ".fna")):
            return ref_path
        if low.endswith(GZIP_FASTA_EXTS):
            return self._gunzip_file(ref_path, "ref_vaf_gunzip")
        raise RuntimeError(f"No FASTA found for VAF in {ref_input}")

    # ------------------------------ helpers: compression (outputs) ------------------------------
    def _pigz_available(self) -> bool:
        """Return True if pigz is enabled and available in PATH."""
        return self.use_pigz.get() and shutil.which("pigz") is not None

    def _gzip_file_fastq_only(self, path: str):
        """
        Gzip a FASTQ file (only if it ends in .fastq or .fq).

        Uses pigz if available (faster on multicore), otherwise Python gzip.
        """
        low = path.lower()
        if not (low.endswith(".fastq") or low.endswith(".fq")):
            return

        try:
            if self._pigz_available():
                args = ["pigz", "-f"]
                if not self.delete_uncompressed_fastq.get():
                    args.append("-k")  # keep original
                args.append(path)
                subprocess.run(args, check=True)
            else:
                gz_path = path + ".gz"
                with open(path, "rb") as src, gzip.open(gz_path, "wb", compresslevel=6) as dst:
                    shutil.copyfileobj(src, dst)
                if self.delete_uncompressed_fastq.get():
                    os.remove(path)

            self.add_log(f"[compress] gzipped {os.path.basename(path)}")
        except Exception as e:
            self.add_log(f"[compress] Failed to gzip {path}: {e}")

    def _compress_fastqs_in(self, d: str):
        """Compress all uncompressed FASTQs in a directory."""
        if not (d and os.path.isdir(d)):
            return
        for pat in ("*.fastq", "*.fq"):
            for f in glob.glob(os.path.join(d, pat)):
                self._gzip_file_fastq_only(f)

    def _maybe_compress_fastqs_now(self, folder: str):
        """Apply compression if the GUI option is enabled."""
        if self.compress_fastq.get():
            self._compress_fastqs_in(folder)

    def _tar_gz_dir(self, dir_path: str, out_name: str | None = None):
        """
        Create a .tar.gz archive of a directory (used for packaging BLAST output, etc.).
        """
        if out_name is None:
            out_name = os.path.basename(os.path.normpath(dir_path)) + ".tar.gz"
        out_path = os.path.join(os.path.dirname(dir_path), out_name)
        try:
            with tarfile.open(out_path, "w:gz") as tar:
                tar.add(dir_path, arcname=os.path.basename(dir_path))
            self.add_log(f"[package] wrote {out_path}")
            return out_path
        except Exception as e:
            self.add_log(f"[package] Failed to tar.gz {dir_path}: {e}")
            return None

    # ------------------------------ process runner ------------------------------
    def add_log(self, message: str):
        """
        Append a log message to the GUI log window in a thread-safe way.

        Tkinter widgets should only be updated on the main thread, so we use .after().
        """
        self.log_text.after(0, lambda: self.log_text.insert(tk.END, message + "\n"))
        self.log_text.after(0, lambda: self.log_text.see(tk.END))

    def run_script_live(self, script_name, *args) -> bool:
        """
        Run a bash script with arguments, streaming output live into the GUI log window.
        """
        cmd = ["bash", script_name] + list(args)
        wrapped = ["stdbuf", "-oL", "-eL"] + cmd
        self.add_log("Running: " + " ".join(shlex.quote(c) for c in wrapped))

        try:
            process = subprocess.Popen(
                wrapped,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True,
            )
        except Exception as e:
            self.add_log(f"❌ Failed to start process {script_name}: {e}")
            return False

        # Stream output line-by-line into the UI
        for line in process.stdout:
            self.add_log(line.rstrip())

        rc = process.wait()
        return rc == 0

    def run_alignment_adaptive(self, script_name, fastq_dir, ref_input, out_dir, threads) -> bool:
        """
        Run alignment with multiple argument styles to support different script conventions.

        Some alignment scripts take:
        - ref directory + '--threads N'
        - ref directory + 'N' as a positional argument
        - or a FASTA file path instead of a directory

        This tries a small set of candidate argument patterns and stops at the first success.
        """
        ref_dir, ref_file = self._ref_dir_and_file(ref_input)

        candidates = [
            ("dir+flag", [fastq_dir, ref_dir, out_dir, "--threads", str(threads)]),
            ("dir+pos",  [fastq_dir, ref_dir, out_dir, str(threads)]),
        ]
        if ref_file:
            candidates += [
                ("file+flag", [fastq_dir, ref_file, out_dir, "--threads", str(threads)]),
                ("file+pos",  [fastq_dir, ref_file, out_dir, str(threads)]),
            ]

        # Remove duplicates while preserving order (so we don't try identical args twice)
        seen, unique = set(), []
        for tag, args in candidates:
            key = tuple(args)
            if key not in seen:
                seen.add(key)
                unique.append((tag, args))

        for tag, args in unique:
            self.add_log(f"[alignment retry] Trying mode: {tag}  args={args}")
            if self.run_script_live(script_name, *args):
                self.add_log(f"[alignment retry] ✅ Succeeded with mode: {tag}")
                return True

        self.add_log("[alignment retry] ❌ All alignment modes failed.")
        return False

    # ------------------------------- UI build -------------------------------
    def build_ui(self):
        """Build the main GUI layout (scrollable steps, buttons, progress bar, logs)."""
        main_frame = tk.Frame(self)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)

        # Scrollable container because the GUI can get quite tall
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas)
        self.scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # One-click section holds all pipeline inputs + step toggles
        self.build_one_click_section(self.scrollable_frame)

        # Separate collapsible panel for launching Streamlit report
        step = CollapsibleStep(
            self.scrollable_frame,
            "interactive_report_streamlit",
            ["Path to app.py", "Port (e.g., 8501)"],
            self.browse_types,
        )
        step.var.set(True)
        step.toggle()
        step.entries["Port (e.g., 8501)"].insert(0, "8501")
        step.pack(fill='x', pady=5)
        self.steps["interactive_report_streamlit"] = step

        # Compression options panel
        fastq_frame = tk.LabelFrame(self.scrollable_frame, text="Compression Options")
        fastq_frame.pack(fill='x', padx=5, pady=6)
        tk.Checkbutton(
            fastq_frame, text="Compress FASTQ outputs to .gz", variable=self.compress_fastq
        ).grid(row=0, column=0, sticky='w', padx=6, pady=4)
        tk.Checkbutton(
            fastq_frame, text="Delete uncompressed FASTQs after compress", variable=self.delete_uncompressed_fastq
        ).grid(row=0, column=1, sticky='w', padx=6)
        tk.Checkbutton(
            fastq_frame, text="Use pigz if available (faster)", variable=self.use_pigz
        ).grid(row=0, column=2, sticky='w', padx=6)
        tk.Checkbutton(
            fastq_frame, text="Also package BLAST output folder (.tar.gz)", variable=self.tar_blast_output
        ).grid(row=1, column=0, columnspan=2, sticky='w', padx=6, pady=(0, 6))

        # Run controls + progress bar
        control_frame = tk.Frame(self)
        control_frame.pack(fill='x', padx=10, pady=10)
        self.progress = ttk.Progressbar(control_frame, orient="horizontal", mode="determinate")
        self.progress.pack(fill='x', pady=5)

        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(pady=5)
        tk.Button(
            btn_frame,
            text="Open Interactive Streamlit Report",
            command=self.open_interactive_report,
            bg="#4CAF50",
            fg="white",
            font=("Arial", 11),
            width=28,
        ).grid(row=0, column=0, padx=6)

        self.run_all_button = tk.Button(
            btn_frame,
            text="Run",
            command=self.start_one_click_pipeline,
            bg="#FF9800",
            fg="white",
            font=("Arial", 11),
            width=12,
        )
        self.run_all_button.grid(row=0, column=1, padx=6)

        tk.Button(
            btn_frame,
            text="Open FastQC Reports",
            command=self.open_fastqc_reports,
            bg="#2196F3",
            fg="white",
            font=("Arial", 11),
            width=20,
        ).grid(row=0, column=2, padx=6)

        tk.Button(
            btn_frame,
            text="Open IGV",
            command=self.launch_igv,
            bg="#9C27B0",
            fg="white",
            font=("Arial", 11),
            width=12,
        ).grid(row=0, column=3, padx=6)

        # Log output area
        self.log_text = scrolledtext.ScrolledText(self, height=12, font=("Arial", 10))
        self.log_text.pack(fill='both', expand=True, padx=10, pady=10)
        self.log_text.insert(tk.END, "Pipeline logs will appear here...\n")

    def build_one_click_section(self, parent):
        """Create the one-click pipeline configuration panel + include-step checkboxes."""
        labels = [
            "Root Output Directory",
            "Initial FASTQ Directory",
            "Human Reference Genome",
            "CMV Reference Genome",
            "BLAST DB Directory",
            "Quality Threshold",
            "Threads/Cores",
            # NEW optional snpEff inputs
            "GFF3 Annotation (for snpEff)",
            "Protein FASTA (optional, snpEff)",
        ]
        step = CollapsibleStep(parent, "run_entire_pipeline (nested folders)", labels, self.browse_types)
        step.var.set(True)
        step.toggle()
        step.entries["Threads/Cores"].insert(0, "8")
        step.entries["Quality Threshold"].insert(0, "20")
        step.pack(fill='x', pady=5)
        self.one_click = step

        toggles_frame = tk.LabelFrame(step.content_frame, text="Include steps")
        toggles_frame.pack(fill='x', padx=5, pady=6)

        # Default “Typical” preset for steps (your choices)
        defaults_true = {
            "split_interleaved": False,
            "filtering": True,
            "index_human": False,
            "align_human": True,
            "unaligned_reads": True,
            "index_cmv": False,
            "align_cmv": True,
            "extract_cmv_metrics": True,
            "blast_summary": False,
            "summarize_blast_hits": False,
            "variant_calling_vaf": True,
            "annotate_vcf_snpeff": True,  # NEW (optional)
        }

        grid_cols = 3
        controls = tk.Frame(toggles_frame)
        controls.grid(row=0, column=0, columnspan=grid_cols, sticky='w', pady=(0, 6))

        # Convenience buttons: tick/untick all, or apply the default preset
        def set_all(val: bool):
            for v in self.one_click_vars.values():
                v.set(val)

        def apply_typical():
            for k, v in defaults_true.items():
                self.one_click_vars[k].set(v)

        tk.Button(controls, text="Select All", command=lambda: set_all(True)).grid(row=0, column=0, padx=4)
        tk.Button(controls, text="Clear All", command=lambda: set_all(False)).grid(row=0, column=1, padx=4)
        tk.Button(controls, text="Typical Preset", command=apply_typical).grid(row=0, column=2, padx=4)

        # Create one checkbox per pipeline step, laid out in a small grid
        for i, name in enumerate(self.one_click_step_order, start=1):
            var = tk.BooleanVar(value=defaults_true.get(name, True))
            cb = tk.Checkbutton(toggles_frame, text=name.replace('_', ' ').title(), variable=var, anchor='w')
            r, c = divmod(i, grid_cols)
            cb.grid(row=r + 1, column=c, sticky='w', padx=6, pady=3)
            self.one_click_vars[name] = var

    # ---------------------------- one-click runner ----------------------------
    def run_one_click_pipeline_worker(self):
        """
        Worker function that executes the selected pipeline steps.

        This runs in a background thread to keep the UI responsive.
        It:
        - reads UI values
        - normalises inputs (archives/gz/single files)
        - creates a timestamped run folder
        - executes each selected step in order
        - updates the progress bar and logs
        """
        self.run_all_button.config(state="disabled")
        try:
            vals = self.one_click.get_values()
            root_dir = vals["Root Output Directory"]
            fastq_start = vals["Initial FASTQ Directory"]
            human_ref = vals["Human Reference Genome"]
            cmv_ref = vals["CMV Reference Genome"]
            blast_db = vals["BLAST DB Directory"]
            qthresh = vals["Quality Threshold"] or "20"
            threads = vals["Threads/Cores"] or "8"
            gff3_path = vals.get("GFF3 Annotation (for snpEff)", "").strip()
            prot_fa = vals.get("Protein FASTA (optional, snpEff)", "").strip()

            # Minimal validation (more checks happen step-by-step)
            if not root_dir:
                raise ValueError("Root Output Directory is required.")
            if not fastq_start:
                raise ValueError("Initial FASTQ Directory is required.")

            # Normalise user-selected inputs into a consistent format for scripts
            fastq_start = self._normalize_fastq_input(fastq_start)
            human_ref = self._normalize_reference_for_alignment(human_ref) if human_ref else ""
            cmv_ref = self._normalize_reference_for_alignment(cmv_ref) if cmv_ref else ""
            blast_db = self._normalize_blast_db(blast_db) if blast_db else ""

            # Each run goes into a timestamp folder (helps traceability & avoids overwriting)
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            run_root = os.path.join(root_dir, timestamp)
            self._ensure_dir(run_root)
            self.add_log(f"[init] Run root = {run_root}")

            current_fastq_dir = self._auto_fastq_dir(fastq_start)
            self._assert_fastqs(current_fastq_dir, "initial")
            self.add_log(f"[init] FASTQ dir = {current_fastq_dir}")

            # These are filled in as the pipeline runs
            human_align_dir = None
            unaligned_fastq_dir = None
            cmv_align_dir = None
            cmv_bam_dir = None
            ref_for_vaf = None  # resolved CMV FASTA used by VAF & snpEff

            # Only run steps the user selected
            selected = [s for s in self.one_click_step_order if self.one_click_vars[s].get()]
            if not selected:
                messagebox.showerror("No steps selected", "Please tick at least one step in 'Include steps'.")
                self.add_log("❌ One-click aborted: no steps selected.")
                return
            self.add_log("One-click selected steps (in order): " + ", ".join(selected))

            # Progress bar is based on number of selected steps
            self.progress["maximum"] = len(selected)
            progressed = 0

            def tick(name):
                """Update progress bar and log a step completion message."""
                nonlocal progressed
                progressed += 1
                self.progress["value"] = progressed
                self.progress.update_idletasks()
                self.add_log(f"✅ Finished: {name}")

            # Run steps in a fixed order (even if user selects a later step)
            for name in self.one_click_step_order:
                if not self.one_click_vars[name].get():
                    continue

                self.add_log(f"▶️ Running: {name}")

                if name == "split_interleaved":
                    # Stage input FASTQs into a working folder to keep layout predictable
                    stage_dir = os.path.join(run_root, "split_input")
                    self._symlink_or_copy_fastqs(current_fastq_dir, stage_dir)
                    ok = self.run_script_live("split_interleaved.sh", stage_dir)
                    if not ok:
                        raise RuntimeError("split_interleaved failed")
                    split_dir = os.path.join(stage_dir, "split_interleaved")
                    current_fastq_dir = split_dir if os.path.isdir(split_dir) else stage_dir
                    self._assert_fastqs(current_fastq_dir, "split_interleaved")
                    tick(name)
                    self._maybe_compress_fastqs_now(current_fastq_dir)

                elif name == "filtering":
                    # Read filtering step (uses your qthresh and threads settings)
                    filt_dir = os.path.join(run_root, "filtering")
                    ok = self.run_script_live("filtering.sh", current_fastq_dir, filt_dir, qthresh, threads)
                    if not ok:
                        raise RuntimeError("filtering failed")
                    current_fastq_dir = os.path.join(filt_dir, "filtered_reads")
                    self._assert_fastqs(current_fastq_dir, "filtering")
                    tick(name)
                    self._maybe_compress_fastqs_now(current_fastq_dir)

                elif name == "index_human":
                    if not human_ref:
                        raise ValueError("Human Reference Genome is required for index_human.")
                    idx_dir = os.path.join(run_root, "filtering", "human_alignment", "index_human")
                    ok = self.run_script_live("index_human.sh", human_ref, idx_dir, threads)
                    if not ok:
                        raise RuntimeError("index_human failed")
                    tick(name)

                elif name == "align_human":
                    if not human_ref:
                        raise ValueError("Human Reference Genome is required for align_human.")
                    human_align_dir = os.path.join(run_root, "filtering", "human_alignment")
                    self._ensure_dir(human_align_dir)
                    ok = self.run_alignment_adaptive(
                        "alignment_bwa-mem.sh",
                        current_fastq_dir,
                        human_ref,
                        human_align_dir,
                        threads,
                    )
                    if not ok:
                        raise RuntimeError("align_human failed")
                    tick(name)

                elif name == "unaligned_reads":
                    # Extract FASTQs that did not align to human (so you can try CMV alignment)
                    if not human_align_dir:
                        human_align_dir = os.path.join(run_root, "filtering", "human_alignment")
                    bam_dir = self._find_bam_dir_from_alignment(human_align_dir)
                    unaligned_base_dir = os.path.join(human_align_dir, "unaligned_reads")
                    ok = self.run_script_live("unaligned_reads.sh", bam_dir, unaligned_base_dir, threads)
                    if not ok:
                        raise RuntimeError("unaligned_reads failed")
                    unaligned_fastq_dir = self._pick_fastq_from_unaligned_out(unaligned_base_dir)
                    self._assert_fastqs(unaligned_fastq_dir, "unaligned_reads")
                    current_fastq_dir = unaligned_fastq_dir
                    tick(name)
                    self._maybe_compress_fastqs_now(unaligned_fastq_dir)

                elif name == "index_cmv":
                    if not cmv_ref:
                        raise ValueError("CMV Reference Genome is required for index_cmv.")
                    idx_dir = os.path.join(
                        run_root,
                        "filtering",
                        "human_alignment",
                        "unaligned_reads",
                        "cmv_alignment",
                        "index_cmv",
                    )
                    ok = self.run_script_live("index_cmv.sh", cmv_ref, idx_dir, threads)
                    if not ok:
                        raise RuntimeError("index_cmv failed")
                    tick(name)

                elif name == "align_cmv":
                    # Align unaligned reads to CMV reference
                    if not cmv_ref:
                        raise ValueError("CMV Reference Genome is required for align_cmv.")
                    cmv_align_dir = os.path.join(
                        run_root,
                        "filtering",
                        "human_alignment",
                        "unaligned_reads",
                        "cmv_alignment",
                    )
                    self._ensure_dir(cmv_align_dir)
                    ok = self.run_alignment_adaptive(
                        "alignment_bwa-mem.sh",
                        current_fastq_dir,
                        cmv_ref,
                        cmv_align_dir,
                        threads,
                    )
                    if not ok:
                        raise RuntimeError("align_cmv failed")
                    cmv_bam_dir = self._find_bam_dir_from_alignment(cmv_align_dir)
                    tick(name)

                elif name == "extract_cmv_metrics":
                    # Generate metrics from CMV BAMs, plus the alignment-based summary CSV used by the report
                    if not cmv_align_dir:
                        raise RuntimeError("extract_cmv_metrics requires CMV alignment BAMs. Enable 'align_cmv' first.")
                    cmv_bam_dir = cmv_bam_dir or self._find_bam_dir_from_alignment(cmv_align_dir)
                    metrics_dir = os.path.join(cmv_align_dir, "output")

                    # Keep your original metrics runner exactly as-is
                    ok = run_bam_metrics(cmv_bam_dir, metrics_dir, log_func=self.add_log)
                    if not ok:
                        raise RuntimeError("extract_cmv_metrics failed")

                    # ✅ NEW: write the CSV that Streamlit expects for HQ mapped + breadth based detection
                    # This overwrites metrics_dir/summary.csv to the correct alignment-metrics format.
                    ok2 = write_cmv_alignment_metrics_csv(
                        bam_dir=cmv_bam_dir,
                        out_dir=metrics_dir,
                        mapq=30,
                        depth_min=1,
                        breadth_method="coverage",
                        log_func=self.add_log,
                    )
                    if not ok2:
                        raise RuntimeError("extract_cmv_metrics failed (could not write alignment summary.csv)")

                    tick(name)

                elif name == "blast_summary":
                    # Optional BLAST-based summary step
                    if not blast_db:
                        raise ValueError("BLAST DB Directory is required for blast_summary.")
                    if not cmv_align_dir:
                        raise RuntimeError(
                            "blast_summary expects CMV alignment folder to route outputs. Enable 'align_cmv' first."
                        )
                    blast_out = os.path.join(cmv_align_dir, "output")
                    self._ensure_dir(blast_out)
                    ok = self.run_script_live("02_blast_cmv_summary.sh", current_fastq_dir, blast_db, blast_out, threads)
                    if not ok:
                        raise RuntimeError("blast_summary failed")
                    tick(name)

                    # Optional packaging of BLAST output directory
                    if self.tar_blast_output.get():
                        bo = os.path.join(blast_out, "blast_output")
                        if os.path.isdir(bo):
                            self._tar_gz_dir(bo)

                elif name == "summarize_blast_hits":
                    # Optional: read BLAST TSVs and produce summary tables
                    if not cmv_align_dir:
                        raise RuntimeError("summarize_blast_hits expects CMV alignment folder. Enable 'align_cmv' first.")
                    blast_tsv_dir = os.path.join(cmv_align_dir, "output", "blast_output")
                    out_dir = os.path.join(cmv_align_dir, "output")
                    self._ensure_dir(out_dir)
                    summarize_blast_hits(blast_tsv_dir, out_dir)
                    tick(name)

                    # Optional packaging of BLAST TSVs
                    if self.tar_blast_output.get():
                        if os.path.isdir(blast_tsv_dir):
                            self._tar_gz_dir(blast_tsv_dir)

                elif name == "variant_calling_vaf":
                    # Variant calling step; also resolves which FASTA is used (important for downstream snpEff)
                    if not cmv_align_dir:
                        raise RuntimeError("variant_calling_vaf requires CMV alignment BAMs. Enable 'align_cmv' first.")
                    ref_for_vaf = self._ref_fasta_file(cmv_ref)
                    cmv_bam_dir = cmv_bam_dir or self._find_bam_dir_from_alignment(cmv_align_dir)
                    vaf_out = os.path.join(cmv_align_dir, "output")
                    self._ensure_dir(vaf_out)
                    ok = self.run_script_live("variant_calling_vaf.sh", cmv_bam_dir, ref_for_vaf, vaf_out, threads)
                    if not ok:
                        raise RuntimeError("variant_calling_vaf failed")
                    tick(name)

                elif name == "annotate_vcf_snpeff":
                    # Optional snpEff build + annotate; only runs if a GFF3 is provided
                    if not cmv_align_dir:
                        raise RuntimeError("annotate_vcf_snpeff expects CMV alignment folder. Enable 'align_cmv' first.")
                    if not gff3_path:
                        self.add_log("[snpEff] Skipping annotation: no GFF3 provided.")
                        tick(name)
                        continue

                    vcf_dir = os.path.join(cmv_align_dir, "output", "vcf")
                    if not os.path.isdir(vcf_dir):
                        raise RuntimeError(f"annotate_vcf_snpeff: VCF folder not found at {vcf_dir}")
                    if not ref_for_vaf:
                        ref_for_vaf = self._ref_fasta_file(cmv_ref)

                    snpeff_work = os.path.join(cmv_align_dir, "output", "snpeff")
                    self._ensure_dir(snpeff_work)

                    ok = self.run_script_live("build_snpeff_db.sh", ref_for_vaf, gff3_path, snpeff_work, prot_fa or "")
                    if not ok:
                        raise RuntimeError("snpEff DB build failed")
                    ok = self.run_script_live("annotate_vcfs_snpeff.sh", vcf_dir, snpeff_work)
                    if not ok:
                        raise RuntimeError("snpEff VCF annotation failed")
                    tick(name)

            self.add_log("🎉 One-click pipeline completed successfully!")
            messagebox.showinfo("Success", "One-click pipeline completed successfully!")

        except Exception as e:
            # Show error in GUI and log it so you have a trace
            messagebox.showerror("One-click Error", str(e))
            self.add_log("One-click Error: " + str(e))
        finally:
            # Always re-enable the Run button and reset progress bar
            self.run_all_button.config(state="normal")
            self.progress["value"] = 0

    def start_one_click_pipeline(self):
        """Start pipeline execution in a daemon thread so the GUI stays responsive."""
        threading.Thread(target=self.run_one_click_pipeline_worker, daemon=True).start()

    # ----------------------------- interactive report launcher -----------------------------
    def open_interactive_report(self):
        """
        Launch the Streamlit report (app.py) on the chosen port.

        If app.py path isn't provided in the UI, it prompts the user to pick it.
        """
        step = self.steps.get("interactive_report_streamlit")
        app_path = None
        port = 8501

        if step:
            vals = step.get_values()
            app_path = vals.get("Path to app.py", "").strip() or None
            p = vals.get("Port (e.g., 8501)", "").strip()
            if p:
                try:
                    port = int(p)
                except ValueError:
                    messagebox.showerror("Error", "Port must be an integer, e.g., 8501")
                    return

        if not app_path:
            app_path = filedialog.askopenfilename(
                title="Select app.py",
                filetypes=[("Python files", "*.py"), ("All files", "*.*")]
            )
            if not app_path:
                return

        ok = self._launch_streamlit(app_path, port)
        if not ok:
            messagebox.showerror("Error", "Failed to launch Streamlit. Check logs.")

    # ----------------------------- misc actions -----------------------------
    def open_fastqc_reports(self):
        """Open FastQC summary.html in a browser (user selects the FastQC output folder)."""
        fastqc_dir = filedialog.askdirectory(title="Select FastQC Output Directory")
        if fastqc_dir:
            index_path = os.path.join(fastqc_dir, "summary.html")
            if os.path.isfile(index_path):
                webbrowser.open_new_tab(f"file://{os.path.abspath(index_path)}")
            else:
                messagebox.showerror("Error", f"Cannot find summary.html in {fastqc_dir}")

    def launch_igv(self):
        """
        Launch IGV for a chosen BAM + reference.

        If the reference is gzipped, it is decompressed to a temp location first.
        """
        try:
            bam_file = filedialog.askopenfilename(title="Select BAM file", filetypes=[("BAM files", "*.bam")])
            if not bam_file:
                return

            ref_file = filedialog.askopenfilename(
                title="Select Reference Genome",
                filetypes=[("FASTA files", "*.fasta *.fa *.fna *.fasta.gz *.fa.gz *.fna.gz"), ("All files", "*.*")]
            )
            if not ref_file:
                return

            if ref_file.lower().endswith(GZIP_FASTA_EXTS):
                ref_file = self._gunzip_file(ref_file, "igv_ref")

            subprocess.Popen(["igv", "-g", ref_file, bam_file])
        except Exception as e:
            messagebox.showerror("Error", f"Failed to launch IGV: {e}")

    def _launch_streamlit(self, app_path: str, port: int) -> bool:
        """
        Start a Streamlit process and pipe its logs into the GUI.

        Returns True if Streamlit is running/already running, otherwise False.
        """
        # Avoid starting multiple Streamlit servers
        if self._streamlit_proc and self._streamlit_proc.poll() is None:
            self.add_log("ℹ️ Streamlit already running.")
            return True

        if not os.path.isfile(app_path):
            self.add_log(f"❌ app.py not found at {app_path}")
            return False

        cmd = ["streamlit", "run", app_path, "--server.headless", "true", "--server.port", str(port)]
        self.add_log("Starting Streamlit: " + " ".join(shlex.quote(c) for c in cmd))

        # CLEAR_CACHE_ON_START can be read by your app.py if you want to clear Streamlit cache on launch
        env = os.environ.copy()
        env["CLEAR_CACHE_ON_START"] = "1"

        try:
            # On Unix: start a new process group so we can terminate the whole group on exit
            preexec = os.setsid if os.name != "nt" else None
            # On Windows: new process group equivalent
            creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0

            self._streamlit_proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                universal_newlines=True,
                env=env,
                preexec_fn=preexec,
                creationflags=creationflags
            )
        except FileNotFoundError:
            self.add_log("❌ Could not find the 'streamlit' executable. Is it installed?")
            return False
        except Exception as e:
            self.add_log(f"❌ Failed to start Streamlit: {e}")
            return False

        # Streamlit logs: pump stdout lines into the GUI
        def _pump_logs():
            try:
                for line in self._streamlit_proc.stdout:
                    self.add_log(line.rstrip())
            except Exception:
                pass

        threading.Thread(target=_pump_logs, daemon=True).start()

        # Try to open the browser tab once the server is likely up
        def _open_when_ready():
            url = f"http://localhost:{port}"
            for _ in range(20):
                time.sleep(0.5)
                try:
                    webbrowser.open_new_tab(url)
                    break
                except Exception:
                    continue
            self.add_log(f"🌐 Interactive report available at: {url}")

        threading.Thread(target=_open_when_ready, daemon=True).start()
        return True

    def on_close(self):
        """
        Clean shutdown:
        - terminate Streamlit if running
        - remove temp folders/files
        - close the GUI
        """
        try:
            if self._streamlit_proc and self._streamlit_proc.poll() is None:
                if os.name == "nt":
                    self._streamlit_proc.terminate()
                else:
                    os.killpg(os.getpgid(self._streamlit_proc.pid), signal.SIGTERM)
        except Exception:
            pass

        self._cleanup_temps()
        self.destroy()


if __name__ == "__main__":
    # Entry point for launching the GUI
    app = PipelineGUI()
    app.mainloop()
