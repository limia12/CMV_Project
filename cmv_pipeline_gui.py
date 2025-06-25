import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
import webbrowser
import threading

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
        self.geometry("850x700")

        self.steps_info = {
            "split_interleaved": ["Input", "Output", "Number of Cores"],
            "filtering": ["Input", "Output", "Quality Threshold", "Number of Cores"],
            "index_human": ["Input", "Output", "Number of Cores"],
            "align_human": ["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"],
            "unaligned_reads": ["Input", "Output", "Number of Cores"],
            "index_cmv": ["Input", "Output", "Number of Cores"],
            "align_cmv": ["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"],
            "blast": ["Forward FASTQ", "Reverse FASTQ", "BLAST DB Directory", "Number of Threads", "Output Directory"]
        }

        self.browse_types = {
            "quality threshold": "text",
            "number of cores": "text",
            "number of threads": "text",
            "reference genome": "dir",
            "forward fastq": "file",
            "reverse fastq": "file",
            "blast db directory": "dir",
            "output directory": "dir",
            "fastq dir": "dir"
        }

        self.steps = {}

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
                  bg="#2196F3", fg="white", font=("Arial", 11), width=15).grid(row=0, column=1, padx=5)

        tk.Button(btn_frame, text="Open IGV", command=self.launch_igv,
                  bg="#9C27B0", fg="white", font=("Arial", 11), width=15).grid(row=0, column=2, padx=5)

        self.log_text = scrolledtext.ScrolledText(self, height=10, font=("Arial", 10))
        self.log_text.pack(fill='both', expand=True, padx=10, pady=10)
        self.log_text.insert(tk.END, "Pipeline logs will appear here...\n")

    def add_log(self, message):
        self.log_text.after(0, lambda: self.log_text.insert(tk.END, message + "\n"))
        self.log_text.after(0, lambda: self.log_text.see(tk.END))

    def run_script_live(self, script_name, *args):
        if "alignment_bwa-mem.sh" in script_name and len(args) == 4:
            args = list(args[:3]) + ["--threads", args[3]]

        cmd = ["bash", script_name] + list(args)
        self.add_log("Running: " + " ".join(cmd))

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            self.add_log(line.rstrip())
        err = process.stderr.read()
        if err:
            self.add_log(err)
        process.wait()
        return process.returncode == 0

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
                    success = self.run_script_live("alignment_bwa-mem.sh", vals["FASTQ Dir"], vals["Reference Genome"], vals["Output"], vals["Number of Cores"])
                elif name == "unaligned_reads":
                    check_fields(["Input", "Output", "Number of Cores"])
                    success = self.run_script_live("unaligned_reads.sh", vals["Input"], vals["Output"], vals["Number of Cores"])
                elif name == "index_cmv":
                    check_fields(["Input", "Output", "Number of Cores"])
                    success = self.run_script_live("index_cmv.sh", vals["Input"], vals["Output"], vals["Number of Cores"])
                elif name == "align_cmv":
                    check_fields(["FASTQ Dir", "Reference Genome", "Output", "Number of Cores"])
                    success = self.run_script_live("alignment_bwa-mem.sh", vals["FASTQ Dir"], vals["Reference Genome"], vals["Output"], vals["Number of Cores"])
                elif name == "blast":
                    check_fields(["Forward FASTQ", "Reverse FASTQ", "BLAST DB Directory", "Number of Threads", "Output Directory"])
                    success = self.run_script_live("BLAST.sh", vals["Forward FASTQ"], vals["Reverse FASTQ"], vals["BLAST DB Directory"], vals["Number of Threads"], vals["Output Directory"])

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
        igv_path = filedialog.askopenfilename(title="Select IGV Executable")
        if igv_path and os.path.isfile(igv_path):
            try:
                subprocess.Popen([igv_path])
            except Exception as e:
                messagebox.showerror("Error", f"Failed to launch IGV: {e}")
        else:
            messagebox.showerror("Error", "Invalid IGV executable path.")

if __name__ == "__main__":
    app = PipelineGUI()
    app.mainloop()