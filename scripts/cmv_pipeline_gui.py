import os
import signal
import subprocess
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
import webbrowser
import glob
import matplotlib
matplotlib.use('Agg')  # headless backend
import shlex
import datetime
import shutil
import tempfile
import tarfile
import zipfile
import gzip

from extract_cmv_bam_metrics import run_bam_metrics  # unchanged
from summarize_blast_hits import summarize_blast_hits  # unchanged

# ---------- archive/fastq constants ----------
ARCHIVE_EXTS = ('.zip', '.tar', '.tar.gz', '.tgz', '.tar.bz2', '.tbz2', '.tar.xz', '.txz')
GZIP_FASTA_EXTS = ('.fa.gz', '.fasta.gz', '.fna.gz')
FASTQ_EXTS = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')


# ------------------------------- UI blocks -------------------------------

class CollapsibleStep(tk.Frame):
    """
    Collapsible block for one pipeline step.
    Shows a checkbox and a panel with inputs and optional Browse buttons.
    """
    def __init__(self, parent, step_name, labels, browse_types=None, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.step_name = step_name
        self.labels = labels
        self.browse_types = browse_types or {}

        # Header checkbox
        self.var = tk.BooleanVar()
        self.checkbox = tk.Checkbutton(
            self,
            text=step_name.replace('_', ' ').title(),
            variable=self.var,
            command=self.toggle,
            font=("Arial", 11, "bold"),
        )
        self.checkbox.pack(anchor='w', pady=2)

        # Content area
        self.content_frame = tk.Frame(self)
        self.entries = {}

        for label in labels:
            row = tk.Frame(self.content_frame)
            lbl = tk.Label(row, text=f"{label}:", width=22, anchor='e')
            ent = tk.Entry(row, width=54)
            browse_type = self.browse_types.get(label.lower(), 'dir')

            if browse_type != 'text':
                btn = tk.Button(row, text="Browse", command=lambda l=label: self.browse_path(l))
                btn.pack(side='right', padx=5)

            lbl.pack(side='left')
            ent.pack(side='left', fill='x', expand=True, padx=5)
            row.pack(fill='x', pady=2)
            self.entries[label] = ent

        # Start collapsed
        self.content_frame.pack_forget()

    def toggle(self):
        if self.var.get():
            self.content_frame.pack(fill='x', padx=25, pady=5)
        else:
            self.content_frame.pack_forget()

    def browse_path(self, label):
        """
        Choose a file/dir based on browse_types.
        For 'dir', if user cancels, we offer an archive picker as a fallback.
        """
        bt = self.browse_types.get(label.lower(), 'dir')
        lbl = label.lower()

        path = ""
        if bt == 'dir':
            path = filedialog.askdirectory(title=f"Select {label}")
            if not path:
                # Archive fallback
                path = filedialog.askopenfilename(
                    title=f"Or select an archived {label}",
                    filetypes=[("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                               ("All files", "*.*")]
                )
        elif bt == 'file':
            if "fastq" in lbl:
                filetypes = [("FASTQ files", "*.fastq *.fq *.fastq.gz *.fq.gz"),
                             ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            elif "fasta" in lbl or "reference" in lbl:
                filetypes = [("FASTA files", "*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz"),
                             ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            elif "bam" in lbl:
                filetypes = [("BAM files", "*.bam"), ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            elif "tsv" in lbl:
                filetypes = [("TSV files", "*.tsv *.tsv.gz"),
                             ("Archives", "*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar.xz *.txz"),
                             ("All files", "*.*")]
            else:
                filetypes = [("All files", "*.*")]
            path = filedialog.askopenfilename(title=f"Select {label}", filetypes=filetypes)
        else:
            return

        if path:
            self.entries[label].delete(0, tk.END)
            self.entries[label].insert(0, path)

    def get_values(self):
        return {label: ent.get().strip() for label, ent in self.entries.items()}


# ------------------------------- Main App -------------------------------

class PipelineGUI(tk.Tk):
    """
    GUI for the CMV pipeline.
    Has per-step panels + a one-click (nested folders) runner.
    """

    def __init__(self):
        super().__init__()
        self.title("CMV Bioinformatics Pipeline GUI")
        self.geometry("900x980")
        self._streamlit_proc = None

        # compression options (for outputs)
        self.compress_fastq = tk.BooleanVar(value=False)              # gzip FASTQ outputs
        self.delete_uncompressed_fastq = tk.BooleanVar(value=False)   # remove *.fastq after gzipping
        self.use_pigz = tk.BooleanVar(value=True)                     # use pigz if available
        self.tar_blast_output = tk.BooleanVar(value=False)            # tar.gz the BLAST output folder

        # Temp registry for unpacked/ staged items
        self._temps = []

        # Per-step panels and fields
        self.steps_info = {
            "split_interleaved": ["Input", "Output", "Number of Cores"],
            "filtering": ["Input", "Output", "Quality Threshold", "Number of Cores"],
            "index_human": ["Input", "Output", "Number of Cores"],
            "align_human": ["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"],
            "unaligned_reads": ["Input", "Output", "Number of Cores"],
            "index_cmv": ["Input", "Output", "Number of Cores"],
            "align_cmv": ["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"],
            "extract_cmv_metrics": ["BAM Directory", "Output Directory"],
            "blast_summary": ["FASTQ Directory", "BLAST DB Directory", "Output Directory", "Threads"],
            "summarize_blast_hits": ["BLAST TSV Directory", "Output Directory"],
            "variant_calling_vaf": ["BAM Directory", "Reference FASTA", "Output Directory", "Threads"],
            "interactive_report_streamlit": ["Path to app.py", "Port (e.g., 8501)"],
        }

        # How each label is browsed
        self.browse_types = {
            "quality threshold": "text",
            "number of cores": "text",

            "reference genome": "dir",
            "reference fasta": "file",

            "forward fastq": "file",
            "reverse fastq": "file",
            "blast db directory": "dir",
            "output directory": "dir",
            "fastq directory": "dir",
            "fastq dir": "dir",
            "bam directory": "dir",
            "blast tsv directory": "dir",
            "csv directory": "dir",
            "path to app.py": "file",
            "port (e.g., 8501)": "text",

            # One-click inputs
            "root output directory": "dir",
            "initial fastq directory": "dir",
            "human reference genome": "dir",
            "cmv reference genome": "dir",
            "threads/cores": "text",
            "path to app.py (optional)": "file",
            "blast db": "dir",
        }

        self.steps = {}
        self.one_click_vars = {}
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
            "interactive_report_streamlit",
        ]

        self.build_ui()
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    # ------------------------------ helpers: temps & archives ------------------------------

    def _register_temp(self, p: str) -> str:
        self._temps.append(p)
        return p

    def _cleanup_temps(self):
        for p in getattr(self, "_temps", []):
            try:
                if os.path.isdir(p):
                    shutil.rmtree(p, ignore_errors=True)
                elif os.path.isfile(p):
                    os.remove(p)
            except Exception:
                pass
        self._temps = []

    def _is_archive(self, p: str) -> bool:
        return bool(p) and os.path.isfile(p) and p.lower().endswith(ARCHIVE_EXTS)

    def _unpack_archive(self, archive_path: str, name_hint: str) -> str:
        """Unpack an archive to a temp dir and return an inner directory to work from."""
        dest = self._register_temp(tempfile.mkdtemp(prefix=f"cmv_{name_hint}_"))
        try:
            shutil.unpack_archive(archive_path, dest)
        except Exception as e:
            raise RuntimeError(f"Failed to unpack {archive_path}: {e}")
        entries = [os.path.join(dest, x) for x in os.listdir(dest)]
        dirs = [d for d in entries if os.path.isdir(d)]
        return dirs[0] if len(dirs) == 1 else dest

    def _gunzip_file(self, gz_path: str, name_hint: str) -> str:
        """Gunzip one file to a temp dir, returning the decompressed path."""
        out_dir = self._register_temp(tempfile.mkdtemp(prefix=f"cmv_{name_hint}_"))
        out_path = os.path.join(out_dir, os.path.basename(gz_path)[:-3])
        with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return out_path

    def _stage_single_file(self, src_file: str, name_hint: str) -> str:
        """Place a single file into a temp folder; returns the folder path (so steps expecting directories still work)."""
        d = self._register_temp(tempfile.mkdtemp(prefix=f"cmv_{name_hint}_"))
        dst = os.path.join(d, os.path.basename(src_file))
        try:
            os.symlink(os.path.abspath(src_file), dst)
        except Exception:
            shutil.copy2(src_file, dst)
        return d

    # ------------------------------ helpers: discovery ------------------------------

    def _ensure_dir(self, p: str):
        os.makedirs(p, exist_ok=True)
        return p

    def _next_subdir(self, parent: str, name: str) -> str:
        path = os.path.join(parent, name)
        os.makedirs(path, exist_ok=True)
        return path

    def _has_fastq(self, path: str) -> bool:
        pats = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
        return any(glob.glob(os.path.join(path, p)) for p in pats)

    def _auto_fastq_dir(self, base: str) -> str:
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
        if not self._has_fastq(d):
            raise RuntimeError(f"{stepname}: no FASTQ files found in {d}.")

    def _ref_dir_and_file(self, ref_path: str):
        if os.path.isdir(ref_path):
            return ref_path, None
        return os.path.dirname(ref_path), ref_path

    def _find_bam_dir_from_alignment(self, align_folder: str) -> str:
        candidate = os.path.join(align_folder, "BAM")
        if os.path.isdir(candidate) and glob.glob(os.path.join(candidate, "*.bam")):
            return candidate
        hits = glob.glob(os.path.join(align_folder, "**", "*.bam"), recursive=True)
        if hits:
            return os.path.dirname(hits[0])
        raise RuntimeError(f"No BAM files found under {align_folder}")

    def _pick_fastq_from_unaligned_out(self, out_dir: str) -> str:
        if self._has_fastq(out_dir):
            return out_dir
        nested = os.path.join(out_dir, "unaligned_reads")
        if os.path.isdir(nested) and self._has_fastq(nested):
            return nested
        return self._auto_fastq_dir(out_dir)

    def _symlink_or_copy_fastqs(self, src_dir: str, dst_dir: str):
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
        """Return a folder containing FASTQs from dir/archive/single file."""
        if os.path.isdir(path):
            return self._auto_fastq_dir(path)
        if self._is_archive(path):
            extracted = self._unpack_archive(path, "fastq")
            return self._auto_fastq_dir(extracted)
        if any(path.lower().endswith(ext) for ext in FASTQ_EXTS):
            return self._stage_single_file(path, "fastq_single")
        raise RuntimeError(f"FASTQ input not found/unsupported: {path}")

    def _normalize_bam_dir(self, path: str) -> str:
        """Return a folder that contains .bam files."""
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
        """BLAST DB must be a directory on disk; extract archive if needed."""
        if os.path.isdir(path):
            return path
        if self._is_archive(path):
            return self._unpack_archive(path, "blastdb")
        raise RuntimeError(f"BLAST DB path must be a folder or an archive: {path}")

    def _normalize_tsv_dir(self, path: str) -> str:
        """Return a folder containing TSVs for BLAST summaries."""
        if os.path.isdir(path):
            return path
        if self._is_archive(path):
            return self._unpack_archive(path, "tsv")
        if path.lower().endswith((".tsv", ".tsv.gz")):
            return self._stage_single_file(path, "tsv_single")
        raise RuntimeError(f"BLAST TSV input not found/unsupported: {path}")

    def _normalize_reference_for_alignment(self, ref: str) -> str:
        """Return either a directory or a FASTA file path (if gz, gunzip) usable by your aligner/indexers."""
        if os.path.isdir(ref):
            return ref
        if self._is_archive(ref):
            extracted = self._unpack_archive(ref, "ref")
            return extracted
        low = ref.lower()
        if low.endswith((".fa", ".fasta", ".fna")):
            return ref
        if low.endswith(GZIP_FASTA_EXTS):
            return self._gunzip_file(ref, "ref_gunzip")
        raise RuntimeError(f"Reference not found/unsupported: {ref}")

    def _ref_fasta_file(self, ref_input: str) -> str:
        """Return a real FASTA path (gunzip/extract if necessary) for callers needing a .fa/.fasta file."""
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
        return self.use_pigz.get() and shutil.which("pigz") is not None

    def _gzip_file_fastq_only(self, path: str):
        """gzip *.fastq / *.fq only; skip *.gz and other formats."""
        low = path.lower()
        if not (low.endswith(".fastq") or low.endswith(".fq")):
            return
        try:
            if self._pigz_available():
                args = ["pigz", "-f"]
                if not self.delete_uncompressed_fastq.get():
                    args.append("-k")
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
        if not (d and os.path.isdir(d)):
            return
        for pat in ("*.fastq", "*.fq"):
            for f in glob.glob(os.path.join(d, pat)):
                self._gzip_file_fastq_only(f)

    def _maybe_compress_fastqs_now(self, folder: str):
        if self.compress_fastq.get():
            self._compress_fastqs_in(folder)

    def _tar_gz_dir(self, dir_path: str, out_name: str | None = None):
        """Create a .tar.gz of dir_path as a sibling file."""
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
        self.log_text.after(0, lambda: self.log_text.insert(tk.END, message + "\n"))
        self.log_text.after(0, lambda: self.log_text.see(tk.END))

    def run_script_live(self, script_name, *args) -> bool:
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

        for line in process.stdout:
            self.add_log(line.rstrip())
        rc = process.wait()
        return rc == 0

    def run_alignment_adaptive(self, script_name, fastq_dir, ref_input, out_dir, threads) -> bool:
        """
        Try several calling conventions until one works.
        Modes:
          dir+flag : script FQ REF_DIR OUT --threads N
          dir+pos  : script FQ REF_DIR OUT N
          file+flag: script FQ REF.fa OUT --threads N
          file+pos : script FQ REF.fa OUT N
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

        seen = set()
        unique = []
        for tag, args in candidates:
            key = tuple(args)
            if key not in seen:
                seen.add(key)
                unique.append((tag, args))

        for tag, args in unique:
            self.add_log(f"[alignment retry] Trying mode: {tag}  args={args}")
            ok = self.run_script_live(script_name, *args)
            if ok:
                self.add_log(f"[alignment retry] ✅ Succeeded with mode: {tag}")
                return True

        self.add_log("[alignment retry] ❌ All alignment modes failed.")
        return False

    # ------------------------------- UI build -------------------------------

    def build_ui(self):
        main_frame = tk.Frame(self)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)

        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas)
        self.scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # One-click section at the top
        self.build_one_click_section(self.scrollable_frame)

        # Per-step collapsible panels
        for step_name, fields in self.steps_info.items():
            step = CollapsibleStep(self.scrollable_frame, step_name, fields, self.browse_types)
            if step_name == "interactive_report_streamlit":
                step.entries["Port (e.g., 8501)"].insert(0, "8501")
            step.pack(fill='x', pady=5)
            self.steps[step_name] = step

        # FASTQ compression options
        fastq_frame = tk.LabelFrame(self.scrollable_frame, text="Compression Options")
        fastq_frame.pack(fill='x', padx=5, pady=6)

        tk.Checkbutton(
            fastq_frame, text="Compress FASTQ outputs to .gz",
            variable=self.compress_fastq
        ).grid(row=0, column=0, sticky='w', padx=6, pady=4)

        tk.Checkbutton(
            fastq_frame, text="Delete uncompressed FASTQs after compress",
            variable=self.delete_uncompressed_fastq
        ).grid(row=0, column=1, sticky='w', padx=6)

        tk.Checkbutton(
            fastq_frame, text="Use pigz if available (faster)",
            variable=self.use_pigz
        ).grid(row=0, column=2, sticky='w', padx=6)

        tk.Checkbutton(
            fastq_frame, text="Also package BLAST output folder (.tar.gz)",
            variable=self.tar_blast_output
        ).grid(row=1, column=0, columnspan=2, sticky='w', padx=6, pady=(0,6))

        # Bottom controls
        control_frame = tk.Frame(self)
        control_frame.pack(fill='x', padx=10, pady=10)

        self.progress = ttk.Progressbar(control_frame, orient="horizontal", mode="determinate")
        self.progress.pack(fill='x', pady=5)

        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(pady=5)

        self.run_button = tk.Button(
            btn_frame, text="Run Pipeline", command=self.start_pipeline,
            bg="#4CAF50", fg="white", font=("Arial", 11), width=15
        )
        self.run_button.grid(row=0, column=0, padx=5)

        self.run_all_button = tk.Button(
            btn_frame, text="Run All (nested folders)", command=self.start_one_click_pipeline,
            bg="#FF9800", fg="white", font=("Arial", 11), width=22
        )
        self.run_all_button.grid(row=0, column=1, padx=5)

        tk.Button(
            btn_frame, text="Open FastQC Reports", command=self.open_fastqc_reports,
            bg="#2196F3", fg="white", font=("Arial", 11), width=18
        ).grid(row=0, column=2, padx=5)

        tk.Button(
            btn_frame, text="Open IGV", command=self.launch_igv,
            bg="#9C27B0", fg="white", font=("Arial", 11), width=15
        ).grid(row=0, column=3, padx=5)

        self.log_text = scrolledtext.ScrolledText(self, height=12, font=("Arial", 10))
        self.log_text.pack(fill='both', expand=True, padx=10, pady=10)
        self.log_text.insert(tk.END, "Pipeline logs will appear here...\n")

    def build_one_click_section(self, parent):
        labels = [
            "Root Output Directory",
            "Initial FASTQ Directory",
            "Human Reference Genome",
            "CMV Reference Genome",
            "BLAST DB Directory",
            "Quality Threshold",
            "Threads/Cores",
            "Path to app.py (optional)",
            "Port (e.g., 8501)",
        ]
        step = CollapsibleStep(parent, "run_entire_pipeline (nested folders)", labels, self.browse_types)
        step.var.set(True)
        step.toggle()
        step.entries["Port (e.g., 8501)"].insert(0, "8501")
        step.entries["Threads/Cores"].insert(0, "8")
        step.entries["Quality Threshold"].insert(0, "20")
        step.pack(fill='x', pady=5)
        self.one_click = step

        # Per-step toggles
        toggles_frame = tk.LabelFrame(step.content_frame, text="Include steps")
        toggles_frame.pack(fill='x', padx=5, pady=6)

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
            "interactive_report_streamlit": False,
        }

        controls = tk.Frame(toggles_frame)
        controls.grid(row=0, column=0, columnspan=3, sticky='w', pady=(0, 6))

        def set_all(val: bool):
            for v in self.one_click_vars.values():
                v.set(val)

        def apply_typical():
            for k, v in defaults_true.items():
                self.one_click_vars[k].set(v)

        tk.Button(controls, text="Select All", command=lambda: set_all(True)).pack(side='left', padx=4)
        tk.Button(controls, text="Clear All", command=lambda: set_all(False)).pack(side='left', padx=4)
        tk.Button(controls, text="Typical Preset", command=apply_typical).pack(side='left', padx=4)

        grid_cols = 3
        for i, name in enumerate(self.one_click_step_order, start=1):
            var = tk.BooleanVar(value=defaults_true.get(name, True))
            cb = tk.Checkbutton(toggles_frame, text=name.replace('_', ' ').title(), variable=var)
            r, c = divmod(i, grid_cols)
            cb.grid(row=r, column=c, sticky='w', padx=6, pady=3)
            self.one_click_vars[name] = var

    # ---------------------------- one-click runner ----------------------------

    def run_one_click_pipeline_worker(self):
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
            app_path = vals["Path to app.py (optional)"]
            port_str = vals["Port (e.g., 8501)"] or "8501"

            if not root_dir:
                raise ValueError("Root Output Directory is required.")
            if not fastq_start:
                raise ValueError("Initial FASTQ Directory is required.")
            try:
                port = int(port_str)
            except ValueError:
                raise ValueError("Port must be an integer, e.g., 8501")

            # Normalize inputs (accept dirs, archives, or single files)
            fastq_start = self._normalize_fastq_input(fastq_start)
            human_ref   = self._normalize_reference_for_alignment(human_ref) if human_ref else ""
            cmv_ref     = self._normalize_reference_for_alignment(cmv_ref)   if cmv_ref else ""
            blast_db    = self._normalize_blast_db(blast_db)                 if blast_db else ""

            # Timestamped run folder
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            run_root = os.path.join(root_dir, timestamp)
            self._ensure_dir(run_root)
            self.add_log(f"[init] Run root = {run_root}")

            # Start from initial FASTQ directory
            current_fastq_dir = self._auto_fastq_dir(fastq_start)
            self._assert_fastqs(current_fastq_dir, "initial")
            self.add_log(f"[init] FASTQ dir = {current_fastq_dir}")

            # Keep track of key locations for downstream steps
            human_align_dir = None
            unaligned_base_dir = None
            unaligned_fastq_dir = None
            cmv_align_dir = None
            cmv_bam_dir = None

            selected = [s for s in self.one_click_step_order if self.one_click_vars[s].get()]
            if not selected:
                messagebox.showerror("No steps selected", "Please tick at least one step in 'Include steps'.")
                self.add_log("❌ One-click aborted: no steps selected.")
                return
            self.add_log("One-click selected steps (in order): " + ", ".join(selected))

            self.progress["maximum"] = len(selected)
            progressed = 0

            def tick(name):
                nonlocal progressed
                progressed += 1
                self.progress["value"] = progressed
                self.progress.update_idletasks()
                self.add_log(f"✅ Finished: {name}")

            # ---- chain through steps with enforced nesting ----
            for name in self.one_click_step_order:
                if not self.one_click_vars[name].get():
                    continue
                self.add_log(f"▶️ Running: {name}")

                if name == "split_interleaved":
                    stage_dir = os.path.join(run_root, "split_input")
                    self._symlink_or_copy_fastqs(current_fastq_dir, stage_dir)
                    ok = self.run_script_live("split_interleaved.sh", stage_dir)
                    if not ok:
                        raise RuntimeError("split_interleaved failed")
                    split_dir = os.path.join(stage_dir, "split_interleaved")
                    if not os.path.isdir(split_dir):
                        raise RuntimeError("split_interleaved did not create expected directory")
                    current_fastq_dir = split_dir
                    self._assert_fastqs(current_fastq_dir, "split_interleaved")
                    tick(name)
                    self._maybe_compress_fastqs_now(current_fastq_dir)

                elif name == "filtering":
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
                    self.add_log(f"[align_human] FASTQ dir = {current_fastq_dir}")
                    self.add_log(f"[align_human] Reference = {human_ref}")
                    ok = self.run_alignment_adaptive("alignment_bwa-mem.sh", current_fastq_dir, human_ref, human_align_dir, threads)
                    if not ok:
                        raise RuntimeError("align_human failed")
                    tick(name)

                elif name == "unaligned_reads":
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
                    idx_dir = os.path.join(run_root, "filtering", "human_alignment", "unaligned_reads", "cmv_alignment", "index_cmv")
                    ok = self.run_script_live("index_cmv.sh", cmv_ref, idx_dir, threads)
                    if not ok:
                        raise RuntimeError("index_cmv failed")
                    tick(name)

                elif name == "align_cmv":
                    if not cmv_ref:
                        raise ValueError("CMV Reference Genome is required for align_cmv.")
                    if not unaligned_fastq_dir:
                        unaligned_fastq_dir = current_fastq_dir
                    cmv_align_dir = os.path.join(run_root, "filtering", "human_alignment", "unaligned_reads", "cmv_alignment")
                    self._ensure_dir(cmv_align_dir)
                    self.add_log(f"[align_cmv] FASTQ dir = {unaligned_fastq_dir}")
                    self.add_log(f"[align_cmv] Reference = {cmv_ref}")
                    ok = self.run_alignment_adaptive("alignment_bwa-mem.sh", unaligned_fastq_dir, cmv_ref, cmv_align_dir, threads)
                    if not ok:
                        raise RuntimeError("align_cmv failed")
                    cmv_bam_dir = self._find_bam_dir_from_alignment(cmv_align_dir)
                    tick(name)

                elif name == "extract_cmv_metrics":
                    if not cmv_align_dir:
                        raise RuntimeError("extract_cmv_metrics requires CMV alignment BAMs. Enable 'align_cmv' first.")
                    cmv_bam_dir = cmv_bam_dir or self._find_bam_dir_from_alignment(cmv_align_dir)
                    metrics_dir = os.path.join(cmv_align_dir, "output")
                    ok = run_bam_metrics(cmv_bam_dir, metrics_dir, log_func=self.add_log)
                    if not ok:
                        raise RuntimeError("extract_cmv_metrics failed")
                    tick(name)

                elif name == "blast_summary":
                    if not blast_db:
                        raise ValueError("BLAST DB Directory is required for blast_summary.")
                    if not unaligned_fastq_dir:
                        unaligned_fastq_dir = current_fastq_dir
                    if not cmv_align_dir:
                        raise RuntimeError("blast_summary expects CMV alignment folder to route outputs. Enable 'align_cmv' first.")
                    blast_out = os.path.join(cmv_align_dir, "output")
                    self._ensure_dir(blast_out)
                    self.add_log(f"[blast_summary] FASTQ dir = {unaligned_fastq_dir}")
                    ok = self.run_script_live("02_blast_cmv_summary.sh", unaligned_fastq_dir, blast_db, blast_out, threads)
                    if not ok:
                        raise RuntimeError("blast_summary failed")
                    tick(name)
                    if self.tar_blast_output.get():
                        bo = os.path.join(blast_out, "blast_output")
                        if os.path.isdir(bo):
                            self._tar_gz_dir(bo)

                elif name == "summarize_blast_hits":
                    if not cmv_align_dir:
                        raise RuntimeError("summarize_blast_hits expects CMV alignment folder to route outputs. Enable 'align_cmv' first.")
                    blast_tsv_dir = os.path.join(cmv_align_dir, "output", "blast_output")
                    out_dir = os.path.join(cmv_align_dir, "output")
                    self._ensure_dir(out_dir)
                    self.add_log(f"[summarize_blast_hits] TSV dir = {blast_tsv_dir}")
                    summarize_blast_hits(blast_tsv_dir, out_dir)
                    tick(name)
                    if self.tar_blast_output.get():
                        if os.path.isdir(blast_tsv_dir):
                            self._tar_gz_dir(blast_tsv_dir)

                elif name == "variant_calling_vaf":
                    if not cmv_align_dir:
                        raise RuntimeError("variant_calling_vaf requires CMV alignment BAMs. Enable 'align_cmv' first.")
                    ref_for_vaf = self._ref_fasta_file(cmv_ref)
                    cmv_bam_dir = cmv_bam_dir or self._find_bam_dir_from_alignment(cmv_align_dir)
                    vaf_out = os.path.join(cmv_align_dir, "output")
                    self._ensure_dir(vaf_out)
                    self.add_log(f"[variant_calling_vaf] BAM dir = {cmv_bam_dir}")
                    self.add_log(f"[variant_calling_vaf] Reference = {ref_for_vaf}")
                    ok = self.run_script_live("variant_calling_vaf.sh", cmv_bam_dir, ref_for_vaf, vaf_out, threads)
                    if not ok:
                        raise RuntimeError("variant_calling_vaf failed")
                    tick(name)

                elif name == "interactive_report_streamlit":
                    if not app_path:
                        raise ValueError("Path to app.py is required for launching the interactive report.")
                    ok = self._launch_streamlit(app_path, port)
                    if not ok:
                        raise RuntimeError("interactive_report_streamlit failed to start")
                    tick(name)

            self.add_log("🎉 One-click pipeline completed successfully!")
            messagebox.showinfo("Success", "One-click pipeline completed successfully!")

        except Exception as e:
            messagebox.showerror("One-click Error", str(e))
            self.add_log("One-click Error: " + str(e))
        finally:
            self.run_all_button.config(state="normal")
            self.progress["value"] = 0

    def start_one_click_pipeline(self):
        threading.Thread(target=self.run_one_click_pipeline_worker, daemon=True).start()

    # --------------------------- per-step runner ---------------------------

    def run_pipeline_worker(self):
        self.run_button.config(state="disabled")
        try:
            selected_steps = [step for step in self.steps.values() if step.var.get()]
            self.progress["maximum"] = len(selected_steps)
            current = 0

            for step in selected_steps:
                name = step.step_name
                vals = step.get_values()

                def check_fields(required):
                    missing = [f for f in required if not vals.get(f)]
                    if missing:
                        raise ValueError(f"Missing inputs for {name}: {', '.join(missing)}")

                self.add_log(f"Running step: {name}")
                success = False

                if name == "split_interleaved":
                    check_fields(["Input", "Output", "Number of Cores"])
                    # Accept dir or archive: normalize to a dir of fastqs
                    input_dir = self._normalize_fastq_input(vals["Input"])
                    success = self.run_script_live("split_interleaved.sh", input_dir)
                    if success:
                        split_dir = os.path.join(input_dir, "split_interleaved")
                        target = split_dir if os.path.isdir(split_dir) else input_dir
                        self._maybe_compress_fastqs_now(target)

                elif name == "filtering":
                    check_fields(["Input", "Output", "Quality Threshold", "Number of Cores"])
                    in_fq = self._normalize_fastq_input(vals["Input"])
                    success = self.run_script_live("filtering.sh", in_fq, vals["Output"], vals["Quality Threshold"], vals["Number of Cores"])
                    if success:
                        out_fastq = os.path.join(vals["Output"], "filtered_reads")
                        target = out_fastq if os.path.isdir(out_fastq) else vals["Output"]
                        self._maybe_compress_fastqs_now(target)

                elif name == "index_human":
                    check_fields(["Input", "Output", "Number of Cores"])
                    ref_in = self._normalize_reference_for_alignment(vals["Input"])
                    success = self.run_script_live("index_human.sh", ref_in, vals["Output"], vals["Number of Cores"])

                elif name == "align_human":
                    check_fields(["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"])
                    fq_dir = self._normalize_fastq_input(vals["FASTQ Dir"])
                    ref    = self._normalize_reference_for_alignment(vals["Reference Genome"])
                    out_dir= vals["Output"]
                    self._assert_fastqs(fq_dir, "align_human")
                    success = self.run_alignment_adaptive("alignment_bwa-mem.sh", fq_dir, ref, out_dir, vals["Number of Cores"])

                elif name == "unaligned_reads":
                    check_fields(["Input", "Output", "Number of Cores"])
                    bam_in = self._normalize_bam_dir(vals["Input"])
                    success = self.run_script_live("unaligned_reads.sh", bam_in, vals["Output"], vals["Number of Cores"])
                    if success:
                        out_fastq = os.path.join(vals["Output"], "unaligned_reads")
                        target = out_fastq if os.path.isdir(out_fastq) else vals["Output"]
                        self._maybe_compress_fastqs_now(target)

                elif name == "index_cmv":
                    check_fields(["Input", "Output", "Number of Cores"])
                    ref_in = self._normalize_reference_for_alignment(vals["Input"])
                    success = self.run_script_live("index_cmv.sh", ref_in, vals["Output"], vals["Number of Cores"])

                elif name == "align_cmv":
                    check_fields(["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"])
                    fq_dir = self._normalize_fastq_input(vals["FASTQ Dir"])
                    ref    = self._normalize_reference_for_alignment(vals["Reference Genome"])
                    out_dir= vals["Output"]
                    self._assert_fastqs(fq_dir, "align_cmv")
                    success = self.run_alignment_adaptive("alignment_bwa-mem.sh", fq_dir, ref, out_dir, vals["Number of Cores"])

                elif name == "extract_cmv_metrics":
                    check_fields(["BAM Directory", "Output Directory"])
                    bam_dir = self._normalize_bam_dir(vals["BAM Directory"])
                    success = run_bam_metrics(bam_dir, vals["Output Directory"], log_func=self.add_log)

                elif name == "blast_summary":
                    check_fields(["FASTQ Directory", "BLAST DB Directory", "Output Directory", "Threads"])
                    fq_dir   = self._normalize_fastq_input(vals["FASTQ Directory"])
                    db_dir   = self._normalize_blast_db(vals["BLAST DB Directory"])
                    out_dir  = vals["Output Directory"]
                    success  = self.run_script_live("02_blast_cmv_summary.sh", fq_dir, db_dir, out_dir, vals["Threads"])
                    if success and self.tar_blast_output.get():
                        bo = os.path.join(out_dir, "blast_output")
                        if os.path.isdir(bo):
                            self._tar_gz_dir(bo)

                elif name == "summarize_blast_hits":
                    check_fields(["BLAST TSV Directory", "Output Directory"])
                    tsv_dir = self._normalize_tsv_dir(vals["BLAST TSV Directory"])
                    summarize_blast_hits(tsv_dir, vals["Output Directory"])
                    success = True
                    if success and self.tar_blast_output.get():
                        if os.path.isdir(tsv_dir):
                            self._tar_gz_dir(tsv_dir)

                elif name == "variant_calling_vaf":
                    check_fields(["BAM Directory", "Reference FASTA", "Output Directory", "Threads"])
                    bam_dir = self._normalize_bam_dir(vals["BAM Directory"])
                    ref_fa  = self._ref_fasta_file(vals["Reference FASTA"])
                    success = self.run_script_live("variant_calling_vaf.sh", bam_dir, ref_fa, vals["Output Directory"], vals["Threads"])

                elif name == "interactive_report_streamlit":
                    check_fields(["Path to app.py", "Port (e.g., 8501)"])
                    app_path = vals["Path to app.py"]
                    try:
                        port = int(vals["Port (e.g., 8501)"])
                    except ValueError:
                        raise ValueError("Port must be an integer, e.g., 8501")
                    success = self._launch_streamlit(app_path, port)

                if not success:
                    self.add_log(f"Step {name} failed. Stopping pipeline.")
                    messagebox.showerror("Step Failed", f"Step '{name}' failed. Check logs for details.")
                    return

                current += 1
                self.progress["value"] = current
                self.progress.update_idletasks()

            self.add_log("Pipeline completed successfully!")
            messagebox.showinfo("Success", "Pipeline completed successfully!")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.add_log("Error: " + str(e))
        finally:
            self.run_button.config(state="normal")
            self.progress["value"] = 0

    # ----------------------------- misc actions -----------------------------

    def start_pipeline(self):
        threading.Thread(target=self.run_pipeline_worker, daemon=True).start()

    def open_fastqc_reports(self):
        fastqc_dir = filedialog.askdirectory(title="Select FastQC Output Directory")
        if fastqc_dir:
            index_path = os.path.join(fastqc_dir, "summary.html")
            if os.path.isfile(index_path):
                webbrowser.open_new_tab(f"file://{os.path.abspath(index_path)}")
            else:
                messagebox.showerror("Error", f"Cannot find summary.html in {fastqc_dir}")

    def launch_igv(self):
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
            # If ref is gz, gunzip to temp first (IGV expects plain fasta or .gz supported? we ensure plain)
            if ref_file.lower().endswith(GZIP_FASTA_EXTS):
                ref_file = self._gunzip_file(ref_file, "igv_ref")
            subprocess.Popen(["igv", "-g", ref_file, bam_file])
        except Exception as e:
            messagebox.showerror("Error", f"Failed to launch IGV: {e}")

    def _launch_streamlit(self, app_path: str, port: int) -> bool:
        if self._streamlit_proc and self._streamlit_proc.poll() is None:
            self.add_log("ℹ️ Streamlit already running.")
            return True
        if not os.path.isfile(app_path):
            self.add_log(f"❌ app.py not found at {app_path}")
            return False

        cmd = ["streamlit", "run", app_path, "--server.headless", "true", "--server.port", str(port)]
        self.add_log("Starting Streamlit: " + " ".join(shlex.quote(c) for c in cmd))

        env = os.environ.copy()
        env["CLEAR_CACHE_ON_START"] = "1"

        try:
            preexec = os.setsid if os.name != "nt" else None
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

        def _pump_logs():
            try:
                for line in self._streamlit_proc.stdout:
                    self.add_log(line.rstrip())
            except Exception:
                pass

        threading.Thread(target=_pump_logs, daemon=True).start()

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
        try:
            if self._streamlit_proc and self._streamlit_proc.poll() is None:
                if os.name == "nt":
                    self._streamlit_proc.terminate()
                else:
                    os.killpg(os.getpgid(self._streamlit_proc.pid), signal.SIGTERM)
        except Exception:
            pass
        # clean temp workspaces
        self._cleanup_temps()
        self.destroy()


if __name__ == "__main__":
    app = PipelineGUI()
    app.mainloop()
