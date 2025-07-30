import os
import subprocess
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
import webbrowser
import pandas as pd
import glob
import matplotlib
matplotlib.use('Agg')  # Use non-interactive Agg backend for Matplotlib (no GUI)

import matplotlib.pyplot as plt
import seaborn as sns
import shlex  # added for safe command logging

from extract_cmv_bam_metrics import run_bam_metrics  # Import the function you need for BAM metrics
from summarize_blast_hits import summarize_blast_hits  # Import the summarize_blast_hits function
import generate_final_report  # Import the final report script

class CollapsibleStep(tk.Frame):
    def __init__(self, parent, step_name, labels, browse_types=None, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.step_name = step_name
        self.labels = labels
        self.browse_types = browse_types or {}

        self.var = tk.BooleanVar()
        self.checkbox = tk.Checkbutton(self, text=step_name.replace('_', ' ').title(),
                                       variable=self.var, command=self.toggle,
                                       font=("Arial", 11, "bold"))
        self.checkbox.pack(anchor='w', pady=2)

        self.content_frame = tk.Frame(self)
        self.entries = {}

        for label in labels:
            row = tk.Frame(self.content_frame)
            lbl = tk.Label(row, text=f"{label}:", width=20, anchor='e')
            ent = tk.Entry(row, width=50)
            browse_type = self.browse_types.get(label.lower(), 'dir')
            if browse_type != 'text':
                btn = tk.Button(row, text="Browse", command=lambda l=label: self.browse_path(l))
                btn.pack(side='right', padx=5)
            lbl.pack(side='left')
            ent.pack(side='left', fill='x', expand=True, padx=5)
            row.pack(fill='x', pady=2)
            self.entries[label] = ent

        self.content_frame.pack_forget()

    def toggle(self):
        if self.var.get():
            self.content_frame.pack(fill='x', padx=25, pady=5)
        else:
            self.content_frame.pack_forget()

    def browse_path(self, label):
        bt = self.browse_types.get(label.lower(), 'dir')
        lbl_lower = label.lower()

        if bt == 'dir':
            path = filedialog.askdirectory(title=f"Select {label}")
        elif bt == 'file':
            path = filedialog.askopenfilename(title=f"Select {label}")
        else:
            if 'reference genome' in lbl_lower:
                path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.*")])
            elif 'fastq' in lbl_lower and ('dir' in lbl_lower or 'directory' in lbl_lower):
                path = filedialog.askdirectory(title=f"Select {label}")
            elif 'fastq' in lbl_lower:
                path = filedialog.askopenfilename(filetypes=[("FASTQ files", "*.fastq *.fq *.fastq.gz *.fq.gz"), ("All files", "*.*")])
            elif 'blast db directory' in lbl_lower:
                path = filedialog.askdirectory(title="Select BLAST DB Directory")
            elif 'csv' in lbl_lower:
                path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
            elif 'bed' in lbl_lower:
                path = filedialog.askopenfilename(filetypes=[("BED files", "*.bed"), ("All files", "*.*")])
            else:
                path = filedialog.askdirectory(title=f"Select {label}")

        if path:
            self.entries[label].delete(0, tk.END)
            self.entries[label].insert(0, path)

    def get_values(self):
        return {label: ent.get().strip() for label, ent in self.entries.items()}


class PipelineGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CMV Bioinformatics Pipeline GUI")
        self.geometry("850x750")

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
            "generate_final_report": ["CSV Directory", "FastQC Directory", "BLAST Directory", "Output Directory"]  # Updated step
        }

        self.browse_types = {
            "quality threshold": "text",
            "number of cores": "text",
            "reference genome": "dir",
            "forward fastq": "file",
            "reverse fastq": "file",
            "blast db directory": "dir",
            "output directory": "dir",
            "fastq dir": "dir",
            "bam directory": "dir",
            "blast tsv directory": "dir",
            "csv directory": "dir",
        }

        self.steps = {}
        self.build_ui()

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

        for step_name, fields in self.steps_info.items():
            step = CollapsibleStep(self.scrollable_frame, step_name, fields, self.browse_types)
            step.pack(fill='x', pady=5)
            self.steps[step_name] = step

        control_frame = tk.Frame(self)
        control_frame.pack(fill='x', padx=10, pady=10)

        self.progress = ttk.Progressbar(control_frame, orient="horizontal", mode="determinate")
        self.progress.pack(fill='x', pady=5)

        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(pady=5)

        self.run_button = tk.Button(btn_frame, text="Run Pipeline", command=self.start_pipeline,
                                    bg="#4CAF50", fg="white", font=("Arial", 11), width=15)
        self.run_button.grid(row=0, column=0, padx=5)

        tk.Button(btn_frame, text="Open FastQC Reports", command=self.open_fastqc_reports,
                  bg="#2196F3", fg="white", font=("Arial", 11), width=18).grid(row=0, column=1, padx=5)

        tk.Button(btn_frame, text="Open IGV", command=self.launch_igv,
                  bg="#9C27B0", fg="white", font=("Arial", 11), width=15).grid(row=0, column=2, padx=5)

        self.log_text = scrolledtext.ScrolledText(self, height=10, font=("Arial", 10))
        self.log_text.pack(fill='both', expand=True, padx=10, pady=10)
        self.log_text.insert(tk.END, "Pipeline logs will appear here...\n")

    def add_log(self, message):
        self.log_text.after(0, lambda: self.log_text.insert(tk.END, message + "\n"))
        self.log_text.after(0, lambda: self.log_text.see(tk.END))

    def run_script_live(self, script_name, *args):
        # UPDATED: merge stderr into stdout and force line-buffering to avoid pipe deadlocks
        cmd = ["bash", script_name] + list(args)
        wrapped = ["stdbuf", "-oL", "-eL"] + cmd

        self.add_log("Running: " + " ".join(shlex.quote(c) for c in wrapped))

        process = subprocess.Popen(
            wrapped,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,   # merge stderr into stdout
            text=True,
            bufsize=1,                  # line-buffered in Python
            universal_newlines=True
        )
        for line in process.stdout:
            self.add_log(line.rstrip())
        rc = process.wait()
        return rc == 0

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
                    success = self.run_script_live("split_interleaved.sh", vals["Input"], vals["Output"], vals["Number of Cores"])
                elif name == "filtering":
                    check_fields(["Input", "Output", "Quality Threshold", "Number of Cores"])
                    success = self.run_script_live("filtering.sh", vals["Input"], vals["Output"], vals["Quality Threshold"], vals["Number of Cores"])
                elif name == "index_human":
                    check_fields(["Input", "Output", "Number of Cores"])
                    success = self.run_script_live("index_human.sh", vals["Input"], vals["Output"], vals["Number of Cores"])
                elif name == "align_human":
                    check_fields(["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"])
                    # UPDATED: pass threads flag correctly
                    success = self.run_script_live(
                        "alignment_bwa-mem.sh",
                        vals["FASTQ Dir"],
                        vals["Reference Genome"],
                        vals["Output"],
                        "--threads", vals["Number of Cores"]
                    )
                elif name == "unaligned_reads":
                    check_fields(["Input", "Output", "Number of Cores"])
                    success = self.run_script_live("unaligned_reads.sh", vals["Input"], vals["Output"], vals["Number of Cores"])
                elif name == "index_cmv":
                    check_fields(["Input", "Output", "Number of Cores"])
                    success = self.run_script_live("index_cmv.sh", vals["Input"], vals["Output"], vals["Number of Cores"])
                elif name == "align_cmv":
                    check_fields(["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"])
                    # UPDATED: pass threads flag correctly
                    success = self.run_script_live(
                        "alignment_bwa-mem.sh",
                        vals["FASTQ Dir"],
                        vals["Reference Genome"],
                        vals["Output"],
                        "--threads", vals["Number of Cores"]
                    )
                elif name == "extract_cmv_metrics":
                    check_fields(["BAM Directory", "Output Directory"])
                    success = run_bam_metrics(vals["BAM Directory"], vals["Output Directory"], log_func=self.add_log)
                elif name == "blast_summary":
                    check_fields(["FASTQ Directory", "BLAST DB Directory", "Output Directory", "Threads"])
                    success = self.run_script_live("02_blast_cmv_summary.sh", vals["FASTQ Directory"], vals["BLAST DB Directory"], vals["Output Directory"], vals["Threads"])
                elif name == "summarize_blast_hits":
                    check_fields(["BLAST TSV Directory", "Output Directory"])  # Updated check
                    blast_tsv_dir = vals["BLAST TSV Directory"]
                    output_dir = vals["Output Directory"]
                    summarize_blast_hits(blast_tsv_dir, output_dir)  # Run the new script function
                    success = True
                elif name == "generate_final_report":
                    check_fields(["CSV Directory", "FastQC Directory", "BLAST Directory", "Output Directory"])
                    csv_dir = vals["CSV Directory"]
                    fastqc_dir = vals["FastQC Directory"]
                    blast_dir = vals["BLAST Directory"]
                    output_dir = vals["Output Directory"]

                    # Find and load the required CSV files from the directory using patterns
                    gene_coverage_csv = glob.glob(os.path.join(csv_dir, "*_sorted_gene.csv"))[0] if glob.glob(os.path.join(csv_dir, "*_sorted_gene.csv")) else None
                    blast_summary_csv = glob.glob(os.path.join(csv_dir, "*_blast_top5.csv"))[0] if glob.glob(os.path.join(csv_dir, "*_blast_top5.csv")) else None
                    bam_metrics_csv = glob.glob(os.path.join(csv_dir, "*_summary.csv"))[0] if glob.glob(os.path.join(csv_dir, "*_summary.csv")) else None
                    per_base_depth_csv = glob.glob(os.path.join(csv_dir, "*_per_base.csv"))[0] if glob.glob(os.path.join(csv_dir, "*_per_base.csv")) else None

                    # Check if all files are found
                    if not all([gene_coverage_csv, blast_summary_csv, bam_metrics_csv, per_base_depth_csv]):
                        raise ValueError("One or more required CSV files are missing in the selected directory.")

                    # Read the CSV files
                    gene_df = pd.read_csv(gene_coverage_csv)
                    blast_df = pd.read_csv(blast_summary_csv)
                    bam_df = pd.read_csv(bam_metrics_csv)
                    per_base_df = pd.read_csv(per_base_depth_csv)

                    # Generate the final report
                    generate_final_report.generate_final_cmv_report(gene_df, blast_df, bam_df, per_base_df, fastqc_dir, blast_dir, output_dir)
                    success = True

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

    def start_pipeline(self):
        threading.Thread(target=self.run_pipeline_worker, daemon=True).start()

    def open_fastqc_reports(self):
        fastqc_dir = filedialog.askdirectory(title="Select FastQC Output Directory")
        if fastqc_dir:
            index_path = os.path.join(fastqc_dir, "summary.html")
            if os.path.isfile(index_path):
                webbrowser.open_new_tab(f"file://{index_path}")
            else:
                messagebox.showerror("Error", f"Cannot find summary.html in {fastqc_dir}")

    def launch_igv(self):
        try:
            bam_file = filedialog.askopenfilename(title="Select BAM file", filetypes=[("BAM files", "*.bam")])
            if not bam_file:
                return
            ref_file = filedialog.askopenfilename(title="Select Reference Genome", filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")])
            if not ref_file:
                return
            subprocess.Popen(["igv", "-g", ref_file, bam_file])
        except Exception as e:
            messagebox.showerror("Error", f"Failed to launch IGV: {e}")


if __name__ == "__main__":
    app = PipelineGUI()
    app.mainloop()
