import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
import webbrowser

class PipelineGUI:
    def __init__(self, root):
        self.root = root
        root.title("CMV Bioinformatics Pipeline GUI")
        root.geometry("800x600")
        root.config(bg="#f0f0f0")

        self.steps = {
            "split_interleaved": tk.BooleanVar(),
            "filtering": tk.BooleanVar(),
            "index_human": tk.BooleanVar(),
            "align_human": tk.BooleanVar(),
            "extract_unaligned": tk.BooleanVar(),
            "index_cmv": tk.BooleanVar(),
            "align_cmv": tk.BooleanVar(),
            "conserved_coverage": tk.BooleanVar(),
        }

        self.paths = {}
        self.build_gui()

    def build_gui(self):
        main_frame = tk.Frame(self.root, bg="#f0f0f0")
        main_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        row = 0
        for step, var in self.steps.items():
            step_frame = tk.Frame(main_frame, bg="#f0f0f0")
            step_frame.grid(row=row, column=0, sticky='w', pady=5)

            tk.Checkbutton(step_frame, text=step.replace('_', ' ').title(),
                           variable=var, bg="#f0f0f0").grid(row=0, column=0, sticky='w', padx=10)

            self.paths[step] = {}
            if step == "conserved_coverage":
                labels = ["Input BAM Dir", "Reference Genome", "Output CSV"]
            elif step in ["align_human", "align_cmv"]:
                labels = ["FASTQ Dir", "Reference Genome", "Output"]
            else:
                labels = ["Input", "Output"]

            for i, label in enumerate(labels):
                entry = tk.Entry(step_frame, width=40)
                entry.grid(row=0, column=1 + i, padx=5)
                button = tk.Button(step_frame, text=f"Browse {label}",
                                   command=lambda s=step, l=label: self.browse_path(s, l))
                button.grid(row=0, column=4 + i, padx=5)
                self.paths[step][label.lower()] = entry
            row += 1

        filter_frame = tk.Frame(main_frame, bg="#f0f0f0")
        filter_frame.grid(row=row, column=0, sticky='w', pady=5)
        tk.Label(filter_frame, text="Enter Quality Threshold:", bg="#f0f0f0").grid(row=0, column=0, padx=10, sticky="w")
        self.paths["filtering"]["quality_threshold"] = tk.Entry(filter_frame, width=40)
        self.paths["filtering"]["quality_threshold"].grid(row=0, column=1, padx=5)

        tk.Label(filter_frame, text="Enter Number of Cores:", bg="#f0f0f0").grid(row=1, column=0, padx=10, sticky="w")
        self.paths["filtering"]["num_cores"] = tk.Entry(filter_frame, width=40)
        self.paths["filtering"]["num_cores"].grid(row=1, column=1, padx=5)
        row += 1

        self.progress = ttk.Progressbar(main_frame, orient="horizontal", mode="determinate", length=600)
        self.progress.grid(row=row, column=0, pady=10)
        row += 1

        button_frame = tk.Frame(main_frame, bg="#f0f0f0")
        button_frame.grid(row=row, column=0, pady=10)
        tk.Button(button_frame, text="Run Pipeline", command=self.run_pipeline,
                  width=20, bg="#4CAF50", fg="white", font=("Arial", 12)).grid(row=0, column=0, padx=5)
        tk.Button(button_frame, text="Open FastQC Reports", command=self.open_fastqc_reports,
                  width=20, bg="#2196F3", fg="white", font=("Arial", 12)).grid(row=0, column=1, padx=5)
        tk.Button(button_frame, text="Open IGV", command=self.launch_igv,
                  width=20, bg="#9C27B0", fg="white", font=("Arial", 12)).grid(row=0, column=2, padx=5)
        row += 1

        self.log_text = scrolledtext.ScrolledText(main_frame, width=90, height=10, wrap=tk.WORD, font=("Arial", 10))
        self.log_text.grid(row=row, column=0, padx=10, pady=10)
        self.log_text.insert(tk.END, "Pipeline logs will appear here...\n")

    def browse_path(self, step, label):
        # For alignment steps, reference genome is a directory, not a file
        if "reference genome" in label.lower():
            # If step is align_human or align_cmv -> pick directory (genome index dir)
            if step in ["align_human", "align_cmv"]:
                path = filedialog.askdirectory(title=f"Select Reference Genome Directory for {step.replace('_', ' ').title()}")
            else:
                # For conserved_coverage or others, pick FASTA file as before
                path = filedialog.askopenfilename(
                    filetypes=[
                        ("FASTA files", "*.fa *.fasta *.fna"),
                        ("All files", "*.*")
                    ]
                )
        elif "csv" in label.lower():
            path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        elif "input" in label.lower():
            path = filedialog.askdirectory()
        else:
            path = filedialog.askdirectory()
        if path:
            self.paths[step][label.lower()].delete(0, tk.END)
            self.paths[step][label.lower()].insert(0, path)

    def run_script(self, script_name, *args):
        try:
            result = subprocess.run(["bash", script_name, *args], capture_output=True, text=True)
            self.log_text.insert(tk.END, result.stdout + "\n")
            self.log_text.yview(tk.END)
            if result.returncode != 0:
                messagebox.showerror("Error", f"Error running {script_name}:\n{result.stderr}")
                return
            return result.stdout
        except Exception as e:
            messagebox.showerror("Error", f"Unexpected error: {e}")

    def run_pipeline(self):
        self.log_text.insert(tk.END, "Running pipeline...\n")
        self.log_text.yview(tk.END)

        total_steps = sum(var.get() for var in self.steps.values())
        current = 0

        shared_cores = self.paths["filtering"]["num_cores"].get()

        def update_progress():
            nonlocal current
            current += 1
            self.progress["value"] = (current / total_steps) * 100
            self.root.update_idletasks()

        if self.steps["split_interleaved"].get():
            self.run_script("split_interleaved_fastq.sh",
                            self.paths["split_interleaved"]["input"].get(),
                            self.paths["split_interleaved"]["output"].get())
            update_progress()

        if self.steps["filtering"].get():
            input_dir = self.paths["filtering"]["input"].get()
            quality_threshold = self.paths["filtering"]["quality_threshold"].get()
            self.run_script("filtering.sh", input_dir, quality_threshold, shared_cores)
            update_progress()

        if self.steps["index_human"].get():
            self.run_script("index_genome.sh",
                            self.paths["index_human"]["input"].get(),
                            self.paths["index_human"]["output"].get())
            update_progress()

        if self.steps["align_human"].get():
            fastq_dir = self.paths["align_human"]["fastq dir"].get()
            reference_genome = self.paths["align_human"]["reference genome"].get()
            output_dir = self.paths["align_human"]["output"].get()
            if shared_cores:
                self.run_script("alignment_bwa-mem.sh", fastq_dir, reference_genome, output_dir, "--threads", shared_cores)
            else:
                self.run_script("alignment_bwa-mem.sh", fastq_dir, reference_genome, output_dir)
            update_progress()

        if self.steps["extract_unaligned"].get():
            self.run_script("unaligned_reads.sh",
                            self.paths["extract_unaligned"]["input"].get(),
                            self.paths["extract_unaligned"]["output"].get())
            update_progress()

        if self.steps["index_cmv"].get():
            self.run_script("index_genome.sh",
                            self.paths["index_cmv"]["input"].get(),
                            self.paths["index_cmv"]["output"].get())
            update_progress()

        if self.steps["align_cmv"].get():
            fastq_dir = self.paths["align_cmv"]["fastq dir"].get()
            reference_genome = self.paths["align_cmv"]["reference genome"].get()
            output_dir = self.paths["align_cmv"]["output"].get()
            if shared_cores:
                self.run_script("alignment_bwa-mem.sh", fastq_dir, reference_genome, output_dir, "--threads", shared_cores)
            else:
                self.run_script("alignment_bwa-mem.sh", fastq_dir, reference_genome, output_dir)
            update_progress()

        if self.steps["conserved_coverage"].get():
            input_bam_dir = self.paths["conserved_coverage"]["input bam dir"].get()
            reference_fasta = self.paths["conserved_coverage"]["reference genome"].get()
            output_csv = self.paths["conserved_coverage"]["output csv"].get()
            if not all([input_bam_dir, reference_fasta, output_csv]):
                messagebox.showwarning("Missing Info", "Please provide all inputs for conserved coverage step.")
            else:
                self.run_script("conserved_coverage.sh", "-i", input_bam_dir, "-r", reference_fasta, "-o", output_csv)
                update_progress()

        messagebox.showinfo("Pipeline Complete", "Pipeline execution is complete.")
        self.log_text.insert(tk.END, "Pipeline execution is complete.\n")
        self.log_text.yview(tk.END)
        self.progress["value"] = 100

    def open_fastqc_reports(self):
        dir_path = filedialog.askdirectory(title="Select FastQC Results Directory")
        if not dir_path:
            return
        for file in os.listdir(dir_path):
            if file.endswith("_fastqc.html"):
                webbrowser.open(f"file://{os.path.abspath(os.path.join(dir_path, file))}")

    def launch_igv(self):
        # Ask user to select reference genome FASTA file
        ref_path = filedialog.askopenfilename(
            title="Select Reference Genome (FASTA)",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.*")]
        )
        if not ref_path:
            return  # User cancelled

        # Ask user to select BAM/VCF files
        file_paths = filedialog.askopenfilenames(
            title="Select BAM/VCF Files to Open in IGV",
            filetypes=[
                ("BAM files", "*.bam"),
                ("VCF files", "*.vcf"),
                ("All files", "*.*")
            ]
        )
        if not file_paths:
            return  # User cancelled

        try:
            # Launch IGV with the reference genome first, then the other files
            subprocess.Popen(["igv", "-g", ref_path, *file_paths])
            self.log_text.insert(tk.END, f"IGV launched with reference: {ref_path} and files: {file_paths}\n")
            self.log_text.yview(tk.END)
        except Exception as e:
            messagebox.showerror("IGV Launch Error", f"Error launching IGV: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = PipelineGUI(root)
    root.mainloop()
