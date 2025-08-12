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

from extract_cmv_bam_metrics import run_bam_metrics  # unchanged
from summarize_blast_hits import summarize_blast_hits  # unchanged
# NOTE: removed: import generate_final_report


class CollapsibleStep(tk.Frame):
    def __init__(self, parent, step_name, labels, browse_types=None, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.step_name = step_name
        self.labels = labels
        self.browse_types = browse_types or {}

        self.var = tk.BooleanVar()
        self.checkbox = tk.Checkbutton(
            self,
            text=step_name.replace('_', ' ').title(),
            variable=self.var,
            command=self.toggle,
            font=("Arial", 11, "bold")
        )
        self.checkbox.pack(anchor='w', pady=2)

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

        self.content_frame.pack_forget()

    def toggle(self):
        if self.var.get():
            self.content_frame.pack(fill='x', padx=25, pady=5)
        else:
            self.content_frame.pack_forget()

    def browse_path(self, label):
        bt = self.browse_types.get(label.lower(), 'dir')
        lbl_lower = label.lower()

        path = ""
        if bt == 'dir':
            path = filedialog.askdirectory(title=f"Select {label}")
        elif bt == 'file':
            path = filedialog.askopenfilename(title=f"Select {label}")
        else:
            # Special cases if you use custom types
            if 'reference genome' in lbl_lower:
                path = filedialog.askopenfilename(filetypes=[
                    ("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.*")
                ], title=f"Select {label}")
            elif 'fastq' in lbl_lower and ('dir' in lbl_lower or 'directory' in lbl_lower):
                path = filedialog.askdirectory(title=f"Select {label}")
            elif 'fastq' in lbl_lower:
                path = filedialog.askopenfilename(filetypes=[
                    ("FASTQ files", "*.fastq *.fq *.fastq.gz *.fq.gz"), ("All files", "*.*")
                ], title=f"Select {label}")
            elif 'blast db directory' in lbl_lower:
                path = filedialog.askdirectory(title="Select BLAST DB Directory")
            elif 'csv' in lbl_lower:
                path = filedialog.askdirectory(title=f"Select {label}")
            elif 'bam stats' in lbl_lower:
                path = filedialog.askdirectory(title="Select BAM Stats Directory (with *_flagstat.txt / *_stats.txt)")
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
        self.geometry("900x780")
        self._streamlit_proc = None

        # Replace "generate_final_report" with "interactive_report_streamlit"
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
            # Final step: launch Streamlit app
            "interactive_report_streamlit": ["Path to app.py", "Port (e.g., 8501)"]
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
            "bam stats directory": "dir",
            "blast tsv directory": "dir",
            "csv directory": "dir",
            "path to app.py": "file",
            "port (e.g., 8501)": "text",
            "reference fasta": "file"
        }

        self.steps = {}
        self.build_ui()
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def build_ui(self):
        main_frame = tk.Frame(self)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)

        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        for step_name, fields in self.steps_info.items():
            step = CollapsibleStep(self.scrollable_frame, step_name, fields, self.browse_types)
            if step_name == "interactive_report_streamlit":
                step.entries["Port (e.g., 8501)"].insert(0, "8501")
            step.pack(fill='x', pady=5)
            self.steps[step_name] = step

        control_frame = tk.Frame(self)
        control_frame.pack(fill='x', padx=10, pady=10)

        self.progress = ttk.Progressbar(control_frame, orient="horizontal", mode="determinate")
        self.progress.pack(fill='x', pady=5)

        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(pady=5)

        self.run_button = tk.Button(
            btn_frame,
            text="Run Pipeline",
            command=self.start_pipeline,
            bg="#4CAF50",
            fg="white",
            font=("Arial", 11),
            width=15
        )
        self.run_button.grid(row=0, column=0, padx=5)

        tk.Button(
            btn_frame,
            text="Open FastQC Reports",
            command=self.open_fastqc_reports,
            bg="#2196F3",
            fg="white",
            font=("Arial", 11),
            width=18
        ).grid(row=0, column=1, padx=5)

        tk.Button(
            btn_frame,
            text="Open IGV",
            command=self.launch_igv,
            bg="#9C27B0",
            fg="white",
            font=("Arial", 11),
            width=15
        ).grid(row=0, column=2, padx=5)

        self.log_text = scrolledtext.ScrolledText(self, height=10, font=("Arial", 10))
        self.log_text.pack(fill='both', expand=True, padx=10, pady=10)
        self.log_text.insert(tk.END, "Pipeline logs will appear here...\n")

    def add_log(self, message):
        self.log_text.after(0, lambda: self.log_text.insert(tk.END, message + "\n"))
        self.log_text.after(0, lambda: self.log_text.see(tk.END))

    def run_script_live(self, script_name, *args):
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
                universal_newlines=True
            )
        except Exception as e:
            self.add_log(f"❌ Failed to start process {script_name}: {e}")
            return False

        for line in process.stdout:
            self.add_log(line.rstrip())
        rc = process.wait()
        return rc == 0

    def _launch_streamlit(self, app_path: str, port: int) -> bool:
        # If already running, do nothing
        if self._streamlit_proc and self._streamlit_proc.poll() is None:
            self.add_log("ℹ️ Streamlit already running.")
            return True

        if not os.path.isfile(app_path):
            self.add_log(f"❌ app.py not found at {app_path}")
            return False

        cmd = [
            "streamlit", "run", app_path,
            "--server.headless", "true",
            "--server.port", str(port),
        ]
        self.add_log("Starting Streamlit: " + " ".join(shlex.quote(c) for c in cmd))

        # Pass env flag so app.py clears Streamlit cache on startup
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
                    success = self.run_script_live(
                        "alignment_bwa-mem.sh",
                        vals["FASTQ Dir"], vals["Reference Genome"], vals["Output"],
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
                    success = self.run_script_live(
                        "alignment_bwa-mem.sh",
                        vals["FASTQ Dir"], vals["Reference Genome"], vals["Output"],
                        "--threads", vals["Number of Cores"]
                    )
                elif name == "extract_cmv_metrics":
                    check_fields(["BAM Directory", "Output Directory"])
                    success = run_bam_metrics(vals["BAM Directory"], vals["Output Directory"], log_func=self.add_log)
                elif name == "blast_summary":
                    check_fields(["FASTQ Directory", "BLAST DB Directory", "Output Directory", "Threads"])
                    success = self.run_script_live("02_blast_cmv_summary.sh", vals["FASTQ Directory"], vals["BLAST DB Directory"], vals["Output Directory"], vals["Threads"])
                elif name == "summarize_blast_hits":
                    check_fields(["BLAST TSV Directory", "Output Directory"])
                    summarize_blast_hits(vals["BLAST TSV Directory"], vals["Output Directory"])
                    success = True
                elif name == "variant_calling_vaf":
                    check_fields(["BAM Directory", "Reference FASTA", "Output Directory", "Threads"])
                    success = self.run_script_live(
                        "variant_calling_vaf.sh",
                        vals["BAM Directory"], vals["Reference FASTA"], vals["Output Directory"], vals["Threads"]
                    )
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
            bam_file = filedialog.askopenfilename(
                title="Select BAM file",
                filetypes=[("BAM files", "*.bam")]
            )
            if not bam_file:
                return
            ref_file = filedialog.askopenfilename(
                title="Select Reference Genome",
                filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
            )
            if not ref_file:
                return
            subprocess.Popen(["igv", "-g", ref_file, bam_file])
        except Exception as e:
            messagebox.showerror("Error", f"Failed to launch IGV: {e}")

    def on_close(self):
        # Cleanly stop Streamlit if it's running
        try:
            if self._streamlit_proc and self._streamlit_proc.poll() is None:
                if os.name == "nt":
                    self._streamlit_proc.terminate()
                else:
                    # only works because we started Streamlit in its own process group
                    os.killpg(os.getpgid(self._streamlit_proc.pid), signal.SIGTERM)
        except Exception:
            pass
        self.destroy()


if __name__ == "__main__":
    app = PipelineGUI()
    app.mainloop()
