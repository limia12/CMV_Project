# CMV Whole Genome Sequencing Pipeline

CMV Whole Genome Sequencing Pipeline is a bioinformatics workflow designed to detect and characterise Cytomegalovirus (CMV) from short-read sequencing data. It combines quality control, alignment, extraction of unmapped reads, BLAST-based confirmation, variant calling, variant annotation, and interactive reporting into one structured workflow.

This pipeline supports both command-line execution and a desktop graphical user interface. It is intended for research and evaluation use.

## Installation

### Prerequisites

- Python version >= 3.10
- Conda installed
- pip installed
- Linux environment recommended
- `environment.yml` available in the project root
- `requirements.txt` available in the project root

### Steps

1. Clone the repository or download the source code:

```bash
git clone <https://github.com/limia12/CMV_Project.git>
cd CMV_Pipeline 
```
2. Create and activate conda environment:
``` bash 
conda env create -f environment.yml
conda activate CMV_Project
```
3. If needed, install any additional python packages from requirements.txt
```bash
pip install -r requirements.txt
```
4. Confirm that the main tools are available:
```bash
bwa
samtools
blastn
lofreq
snpEff
streamlit
```
5. Prepare the reference genome:
```bash 
bash index_genome.sh reference.fa
```
6. Build the snpEff database:
```bash 
bash build_snpeff_db.sh reference.fa annotation.gff snpeff_db/
```

7. Create the BLAST database:
```bash
makeblastdb -in viral_sequences.fna -dbtype nucl -out viral_nt_db
```
## Usage
1. Quality Control and Filtering
    Run filtering.sh to trim paired-end FASTQ files and generate FastQC reports:
```bash 
bash filtering.sh input_fastq/ qc_output/ 20 4
```
This step will:

    Trim reads using Trim Galore
    Generate FastQC reports
    Save filtered reads into qc_output/filtered_reads/
    Save FastQC results into qc_output/fastqc_results/

    Arguments:

    input_fastq/: directory or archive containing input FASTQ files
    qc_output/: output directory
    20: quality threshold
    4: number of cores
    Project Structure

2. Alignment
    Run alignment_bwa-mem.sh to align paired-end reads to the CMV reference genome:
```bash 
bash alignment_bwa-mem.sh qc_output/filtered_reads reference.fa output/ --threads 4
```
    This step will:

    Align reads to the CMV reference genome
    Produce sorted BAM files
    Index BAM files
    Generate BAM statistics

    Outputs are written to:

    output/BAM/
    output/BAM/BAM_Stats/

3. Extract Unaligned Reads
    Run unaligned_reads.sh to extract read pairs where both mates are unmapped:
```bash
bash unaligned_reads.sh output/BAM output/ 4
```
    This step will:

    Process BAM files in the BAM directory
    Extract read pairs where both reads are unmapped
    Write FASTQ outputs into output/unaligned_reads/

4. Split Interleaved FASTQ Files
    If your input FASTQ files are interleaved, run:
```bash
bash split_interleaved_fastq.sh input_fastq/
```
    This step will:

    Split each interleaved FASTQ file into separate R1 and R2 files
    Save the outputs into input_fastq/split_interleaved/
5. BLAST Analysis

    Run 02_blast_cmv_summary.sh to convert FASTQ files to FASTA, run BLAST, and summarise the results:
``` bash
bash 02_blast_cmv_summary.sh input_fastq/ blast_db/ blast_results/ 4
```
    This step will:

    Convert FASTQ files to FASTA
    Run BLAST against the local viral database
    Produce per-sample BLAST TSV files
    Produce summary CSV files

    Outputs are written to:

    blast_results/fasta/
    blast_results/blast_output/
    blast_results/
6. Variant Calling

    Run variant_calling_vaf.sh to call variants and calculate variant allele frequencies:

```bash
bash variant_calling_vaf.sh output/BAM reference.fa variant_results/ 4
```

    This step will:

    Run LoFreq on sorted BAM files
    Create VCF files
    Create VAF tables containing depth and allele frequency values

    Outputs are written to:

    variant_results/vcf/
    variant_results/vaf_tables/

7. Variant Annotation
    Run annotate_vcfs_snpeff.sh to annotate VCF files using snpEff:

```bash 
bash annotate_vcfs_snpeff.sh variant_results/vcf snpeff_db/ 4
```

This step will:

Read VCF files from the VCF directory
Annotate variants using the snpEff database
Write annotated VCF outputs

8. Generate BAM Coverage Metrics
    Run the BAM metrics script to generate summary, gene-level, and per-base coverage outputs:

```bash
python extract_cmv_bam_metrics.py
```

    This step generates:

    <sample>_summary.csv
    <sample>_gene.csv
    <sample>_per_base.csv

    These outputs include:

    Total read count
    Average genome depth
    Gene-level average depth
    Breadth of coverage
    Per-base depth values

9. Summarise BLAST Hits
    Run the Python BLAST summary script:

```bash
python summarize_blast_hits.py --input_dir blast_results/blast_output --output_dir blast_results
```

This step generates:

    <sample>_blast_top5.csv
    <sample>_blast_heatmap.csv
    blast_heatmap.csv

10. Run the Interactive Report
    Launch the Streamlit report:

```bash
streamlit run app.py
```

    The report can display:

    Gene-level depth summaries
    Per-base depth tables
    Alignment summary metrics
    BLAST summaries
    Raw BLAST hits
    Variant allele frequency tables
    FastQC reports
    Annotated VCF results
11. Run the Desktop GUI
    Launch the graphical user interface:

```bash 
python cmv_pipeline_gui.py
```

    The GUI allows you to:

    Select input FASTQ folders or archives
    Select reference genome and BLAST database
    Run selected pipeline steps
    View logs in real time
    Generate alignment metrics CSV files
    Launch the Streamlit report
    Open FastQC output
    Launch IGV for review
    Example Workflow
    bash filtering.sh input_fastq/ qc_output/ 20 4
    bash alignment_bwa-mem.sh qc_output/filtered_reads reference.fa output/ --threads 4
    bash unaligned_reads.sh output/BAM output/ 4
    bash 02_blast_cmv_summary.sh input_fastq/ blast_db/ blast_results/ 4
    bash variant_calling_vaf.sh output/BAM reference.fa variant_results/ 4
    bash annotate_vcfs_snpeff.sh variant_results/vcf snpeff_db/ 4
    streamlit run app.py

CMV_Pipeline/
├── scripts/
│   ├── filtering.sh
│   ├── alignment_bwa-mem.sh
│   ├── unaligned_reads.sh
│   ├── 02_blast_cmv_summary.sh
│   ├── variant_calling_vaf.sh
│   ├── annotate_vcfs_snpeff.sh
│   ├── build_snpeff_db.sh
│   ├── index_genome.sh
│   ├── split_interleaved_fastq.sh
│   ├── extract_cmv_bam_metrics.py
│   ├── summarize_blast_hits.py
│   ├── app.py
├── outputs/
│   ├── BAM/
│   ├── blast_output/
│   ├── vcf/
│   ├── vaf_tables/
├── README.md

## Description of Key Components
    filtering.sh: trims paired-end reads and runs FastQC
    alignment_bwa-mem.sh: aligns reads to the reference genome and generates BAM outputs
    unaligned_reads.sh: extracts read pairs where both mates are unmapped
    split_interleaved_fastq.sh: splits interleaved FASTQ files into separate R1 and R2 files
    02_blast_cmv_summary.sh: runs BLAST on FASTQ-derived FASTA files and produces summary outputs
    summarize_blast_hits.py: summarises BLAST TSV outputs into CSV files
    variant_calling_vaf.sh: calls variants with LoFreq and produces VAF tables
    build_snpeff_db.sh: builds a local snpEff database from a reference FASTA and GFF3 file
    annotate_vcfs_snpeff.sh: annotates VCF files using snpEff
    extract_cmv_bam_metrics.py: calculates BAM summary, gene-level, and per-base coverage metrics
    app.py: interactive Streamlit report for reviewing pipeline outputs
    cmv_pipeline_gui.py: desktop GUI runner for the full pipeline
    environment.yml: conda environment definition for the project
    requirements.txt: additional Python package requirements
    Logs and Outputs

    Outputs are written to the directories specified when running each script. Typical outputs include:

    Filtered FASTQ files
    FastQC reports
    BAM files and BAM statistics
    Unaligned FASTQ files
    BLAST TSV and CSV summaries
    VCF files
    VAF tables
    Annotated VCF files
    Coverage summary CSV files
    Interactive Streamlit visualisations
    Common Issues
    Tool not found: ensure the conda environment is activated and required software is installed
    Reference indexing errors: ensure the reference FASTA exists and indexing completed successfully
    BLAST database not found: ensure makeblastdb was run and viral_nt_db.* files exist
    No VCF files produced: check that BAM inputs are sorted and aligned correctly
    Empty BLAST outputs: check input FASTQ quality and database validity
    GUI launch issues: ensure tkinter and Python dependencies are installed in the active environment
    Streamlit report not loading data: confirm that the expected output files were generated