README.md
# CMV Whole Genome Sequencing Pipeline

The CMV Whole Genome Sequencing Pipeline is designed to detect and characterise Cytomegalovirus (CMV) from short-read sequencing data. It integrates multiple bioinformatics steps into a structured workflow to generate interpretable outputs including viral detection metrics, BLAST-based confirmation, variant calling, and annotation.

This pipeline is intended for research and evaluation purposes, particularly in transplant and virology settings, where detection sensitivity and interpretation of viral signal are important.

---

## Installation

### Prerequisites

- Python ≥ 3.8  
- Conda (recommended)  
- pip installed  

### Required Tools

- bwa  
- samtools  
- fastqc  
- trim_galore  
- blastn (NCBI BLAST+)  
- lofreq  
- snpEff  

---

### Steps

#### 1. Clone the repository

```bash
git clone <your_repository_url>
cd CMV_Pipeline
2. Create and activate environment
conda create -n cmv_pipeline python=3.10
conda activate cmv_pipeline
3. Install dependencies
conda install -c bioconda bwa samtools fastqc trim-galore blast lofreq snpeff
conda install -c conda-forge pandas numpy matplotlib streamlit plotly pysam
4. Prepare reference genome
bash index_genome.sh reference.fa
5. Build snpEff database
bash build_snpeff_db.sh reference.fa annotation.gff snpeff_db/
6. Create BLAST database
makeblastdb -in viral_sequences.fna -dbtype nucl -out viral_nt_db
Usage
1. Quality Control and Filtering
bash filtering.sh input_fastq/ qc_output/ 20 4

This step:

Trims reads using Trim Galore
Generates FastQC reports
2. Alignment
bash alignment_bwa-mem.sh qc_output/ reference.fa output/ --threads 4

This step:

Aligns reads to CMV reference
Produces sorted BAM files
Generates alignment statistics
3. Extract Unaligned Reads
bash unaligned_reads.sh output/BAM results/

This extracts reads where both mates are unmapped for further analysis.

4. BLAST Analysis
bash 02_blast_cmv_summary.sh input_fastq/ blast_db/ results/ 4

This step:

Converts FASTQ to FASTA
Runs BLAST against viral database
Produces summary CSV files
5. Variant Calling
bash variant_calling_vaf.sh output/BAM reference.fa results/

This step:

Calls variants using LoFreq
Calculates allele frequencies
6. Variant Annotation
bash annotate_vcfs_snpeff.sh results/vcf snpeff_db/

This step:

Annotates variants
Identifies potential resistance mutations
7. Visualisation
streamlit run app.py

This launches an interactive report to explore results.

Pipeline Outputs

The pipeline generates several output files:

File	Description
*_summary.csv	Alignment-based CMV metrics
*_gene.csv	Gene-level coverage
*_per_base.csv	Per-base depth values
*_blast.tsv	Raw BLAST hits
*_blast_top5.csv	Top species matches
*_vaf.tsv	Variant allele frequencies
*.annot.vcf	Annotated variants
Detection Approach

CMV detection is based on three complementary metrics:

Reads per million (RPM) – reflects viral abundance
Genome breadth – reflects how much of the genome is covered
BLAST identity – confirms sequence specificity

These metrics are combined to improve sensitivity and specificity and reduce false positives.

Project Structure
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
Description of Key Components
filtering.sh
Performs read trimming and quality control.
alignment_bwa-mem.sh
Aligns reads to the CMV reference genome and produces BAM files.
02_blast_cmv_summary.sh
Runs BLAST against a viral database and summarises species matches.
variant_calling_vaf.sh
Calls variants and calculates allele frequencies.
annotate_vcfs_snpeff.sh
Adds functional annotations to variants.
app.py
Streamlit application for visualising pipeline outputs.
Logs and Outputs
Intermediate and final outputs are written to the specified output directory
FastQC reports are generated for quality assessment
CSV files summarise key results for downstream interpretation
Common Issues
Missing tools
Ensure all dependencies are installed and available in PATH
BLAST performance issues
Increase number of threads
Alignment errors
Ensure reference genome is properly indexed
Empty outputs
Check input FASTQ quality and read counts
Limitations
Developed using synthetic datasets
Requires validation on real clinical samples
Dependent on quality of reference genome and database