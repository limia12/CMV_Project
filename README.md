# CMV_Pipeline_Project

This project helps you analyze sequencing data for **Cytomegalovirus (CMV)**.  
It guides you through the different steps of the pipeline using a simple window you can click through (no coding needed).  
At the end, it creates an **interactive report** you can open in your web browser to explore your results.

---

## 1. Installation

### What you need first
- **Conda** installed (Anaconda or Miniconda)  
- **Python 3.11 or higher**  
- A computer with a desktop (so the window can open)  
- Some bioinformatics tools available (usually through Bioconda):  
  - `bwa`, `samtools`, `bcftools`, `seqtk`, `blast`, `fastqc`, `trim-galore`, `lofreq`, `parallel`, `picard`

### Set up the project

1. Download the project (clone or unzip):
   ```bash
   git clone https://github.com/<your-username>/CMV_Pipeline_Project.git
   cd CMV_Pipeline_Project

2. Create and activate the environment:
    ```bash
    conda env create -f environment.yml
    conda activate cmv_pipeline

3. Install the Python packages:

    pip install -r requirements.txt

4. Running the Application

### To start the program, type:

- python scripts/cmv_pipeline_gui.py

- A window will open.

- This is where you choose what parts of the pipeline you want to run.

- Each step can be turned on or off with a checkbox.

- For most steps you’ll just need to browse for the input folder and output folder, and sometimes pick the number of cores (threads) to use.

- When you’re ready, click Run Pipeline.

- The program will run the steps you selected and show progress in the window.

- At the end, it will open your interactive report in a browser.

## 5. What Each Step Does

Here’s a short explanation of the options you’ll see in the window:

### Split interleaved FASTQs

- Use this if your FASTQ files have both read pairs in one file.

- Input: folder with interleaved FASTQs

- Output: folder where R1 and R2 files will be saved

### Filtering (quality trimming + FastQC)

- Cleans and trims your raw reads

- Input: folder with FASTQs

- Output: will create filtered_reads/ and fastqc_results/

### Index Human Genome (optional)

- Prepares the human genome reference for alignment

- Input: folder with the reference FASTA

- Output: folder to save the index

### Align to Human Genome (optional)

- Aligns your reads to the human genome (to filter human reads)

- Input: FASTQ folder + indexed reference

- Output: BAM files + BAM statistics

### Unaligned Reads (optional)

- Extracts the reads that did not align to the human genome

### Index CMV Genome

- Same as human indexing, but for CMV

### Align to CMV Genome

- Aligns your reads to the CMV reference

- Output: BAM files and alignment stats

### BLAST

- Runs BLAST on the reads to check for CMV and other viruses

### Variant Calling (LoFreq)

- Finds small variants (SNPs/indels) and makes a table of variant allele frequencies

### Generate Report

- Collects all the outputs and opens the interactive report in your browser

---

## 6. The Interactive Report

When the pipeline finishes, your browser will open with the **CMV Interactive Report**.  
This report has several sections (you can expand and collapse them):

### Alignment Statistics
- How many reads mapped vs unmapped  
- Percentage properly paired  
- Average insert size  

### Coverage per Gene
- Bar chart showing average depth across key CMV genes  

### BLAST Results
- Top species detected in your data  
- Heatmaps showing BLAST identity across the genome and within conserved regions  

### Per-base Coverage
- Heatmaps for depth of coverage within specific CMV regions  

### Genome-wide Coverage
- Line plot of depth across the whole CMV genome  

### FastQC Reports
- Summary tables for read quality checks  
- Option to preview and download the full FastQC reports  

### Variants (LoFreq VAF)
- Table of variants with allele frequencies  
- Plots showing allele frequency distribution and positions  
- A simple flag to suggest if there may be mixed strains or co-infection  


7. Project Structure

CMV_Pipeline_Project/
├── scripts/
│   ├── cmv_pipeline_gui.py          # Main window application
│   ├── app.py                       # Interactive report
│   ├── alignment_bwa-mem.sh         # Alignment script
│   ├── filtering.sh                 # Trimming and FastQC
│   ├── index_genome.sh              # Indexing
│   ├── split_interleaved_fastq.sh   # Splitting interleaved reads
│   ├── unaligned_reads.sh           # Extract unmapped reads
│   ├── variant_calling_vaf.sh       # Variant calling
│   ├── 02_blast_cmv_summary.sh      # BLAST
│   ├── extract_cmv_bam_metrics.py   # Coverage metrics
│   └── summarize_blast_hits.py      # BLAST summaries
├── fastq/                           # Example FASTQs
├── output/                          # Outputs: BAMs, CSVs, VCFs, plots
├── References/                      # Reference genomes
├── Simulated_Reads/                 # Simulated data
├── TEST/                            # Small test datasets
├── environment.yml                  # Conda environment
└── requirements.txt                 # Python dependencies

8. Tips if You Get Stuck

    If nothing appears in the report, check that your CSV and BLAST folders are set correctly.

    If FastQC previews look odd, use the download button to open them locally.

    If the report looks empty, try lowering the “BLAST min % identity” in the sidebar.

    To reset things, close the browser tab, go back to your terminal, and restart the app.