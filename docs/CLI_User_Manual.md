# User Guide for CMV Whole Genome Sequencing Pipeline CLI

## Overview

The CMV Whole Genome Sequencing Pipeline is a command-line pipeline designed to detect and characterise Cytomegalovirus (CMV) from short-read sequencing data.

The CLI tools allow users to:

- Trim and quality check FASTQ files
- Align reads to a CMV reference genome
- Extract unmapped reads
- Run BLAST to confirm viral sequence identity
- Call variants and calculate variant allele frequencies
- Annotate variants using snpEff
- Generate files that can be reviewed in the interactive report

This guide explains how to run each command in a simple step-by-step way.

---

## General Usage

Each step of the pipeline is run as a separate command.

Most users will run the steps in this order:

1. Filter and quality check reads
2. Align reads to the CMV reference
3. Extract unmapped reads if needed
4. Run BLAST analysis
5. Call variants
6. Annotate variants
7. Review results in the report

---

## 1. Filter and Quality Check Reads

Use `filtering.sh` to trim paired-end FASTQ files and run FastQC.

### Command

```bash
bash filtering.sh <input_dir_or_archive> <output_dir> <quality_threshold> <num_cores>
```

Example
```bash
bash filtering.sh input_fastq/ qc_output/ 20 4
```
What this script does

This script will:

Read FASTQ files from a folder or supported archive
Trim reads using Trim Galore
Run FastQC on the reads
Save filtered reads
Save FastQC reports
Input arguments
input_dir_or_archive: folder or archive containing FASTQ files
output_dir: folder where output files will be saved
quality_threshold: trimming quality threshold
num_cores: number of CPU cores to use
Output

Typical outputs include:

```bash
filtered_reads/
fastqc_results/
```

2. Align Reads to the CMV Reference Genome

Use alignment_bwa-mem.sh to align paired-end FASTQ files to the CMV reference genome.

Command
```bash
bash alignment_bwa-mem.sh <fastq_dir_or_archive> <genome_index_dir_or_fasta> <output_root> [--threads N]
```

Example
```bash
bash alignment_bwa-mem.sh qc_output/filtered_reads reference.fa output/ --threads 4
```
What this script does

This script will:

Read FASTQ files from a folder or archive
Align paired-end reads using BWA-MEM
Sort BAM files
Index BAM files
Write BAM statistics
Input arguments
fastq_dir_or_archive: folder or archive containing FASTQ files
genome_index_dir_or_fasta: indexed genome directory or reference FASTA
output_root: main output directory
--threads N: optional number of threads
Output

Typical outputs include:

```bash
output/BAM/
output/BAM/BAM_Stats/
```

3. Extract Unmapped Reads

Use unaligned_reads.sh to extract read pairs where both reads are unmapped.

Command
```bash
bash unaligned_reads.sh <bam_dir> <out_base> [threads]
```
Example
```bash
bash unaligned_reads.sh output/BAM output/ 4
```
What this script does

This script will:

Read BAM files from the BAM directory
Find pairs where both mates are unmapped
Convert these reads back into paired FASTQ files
Input arguments
bam_dir: directory containing BAM files
out_base: base output directory
threads: optional number of threads
Output

Typical outputs include:

```bash 
output/unaligned_reads/<sample>_R1.fastq
output/unaligned_reads/<sample>_R2.fastq
``` 
4. Split Interleaved FASTQ Files

Use split_interleaved_fastq.sh if your FASTQ files are interleaved rather than already split into R1 and R2.

Command
```bash
bash split_interleaved_fastq.sh <input_dir>
```
Example
```bash
bash split_interleaved_fastq.sh input_fastq/
```
What this script does

This script will:

Read interleaved FASTQ files
Split them into two files:
read 1
read 2
Input arguments
input_dir: directory containing interleaved FASTQ files
Output

Typical outputs include:

```bash
split_interleaved/<sample>_R1.fastq
split_interleaved/<sample>_R2.fastq
```
5. Run BLAST Analysis

Use 02_blast_cmv_summary.sh to convert FASTQ files to FASTA, run BLAST, and summarise the results.

Command
```bash
bash 02_blast_cmv_summary.sh <FASTQ_DIR_OR_ARCHIVE> <BLAST_DB_DIR> <OUT_DIR> <THREADS>
```
Example
```bash
bash 02_blast_cmv_summary.sh input_fastq/ blast_db/ blast_results/ 4
```
What this script does

This script will:

Read FASTQ files from a folder or archive
Convert FASTQ files to FASTA
Run BLAST against the local viral database
Save per-sample BLAST hits as TSV files
Run the summary script to create CSV summary files
Input arguments
FASTQ_DIR_OR_ARCHIVE: FASTQ folder or archive
BLAST_DB_DIR: folder containing the BLAST database
OUT_DIR: output directory
THREADS: number of BLAST threads
Output

Typical outputs include:

```bash
OUT_DIR/fasta/
OUT_DIR/blast_output/
*_blast.tsv
*_blast_top5.csv
blast_heatmap.csv
```
6. Call Variants

Use variant_calling_vaf.sh to call variants with LoFreq and calculate allele frequencies.

Command
```bash
bash variant_calling_vaf.sh <bam_dir> <ref_fasta(.gz)> <output_dir> [threads]
```
Example
```bash
bash variant_calling_vaf.sh output/BAM reference.fa variant_results/ 4
```

What this script does

This script will:

Read sorted BAM files
Run LoFreq variant calling
Write VCF files
Create TSV files containing:
total depth
alternate depth
allele frequency
Input arguments
bam_dir: directory containing sorted BAM files
ref_fasta(.gz): reference FASTA file
output_dir: output directory
threads: optional number of threads
Output

Typical outputs include:

```bash
variant_results/vcf/
variant_results/vaf_tables/
```

7. Annotate Variants

Use annotate_vcfs_snpeff.sh to annotate VCF files using snpEff.

Command
```bash
bash annotate_vcfs_snpeff.sh <vcf_dir> <snpeff_work_dir> [threads]
```

Example
```bash
bash annotate_vcfs_snpeff.sh variant_results/vcf snpeff_db/ 4
```
What this script does

This script will:

Read VCF files from a folder
Add snpEff annotations
Write annotated VCF outputs
Input arguments
vcf_dir: folder containing VCF files
snpeff_work_dir: folder containing the snpEff database and config
threads: optional compression threads
Output

Typical outputs include:
```bash
*.annot.vcf
*.annot.vcf.gz
```
8. Generate Coverage Metrics

Use extract_cmv_bam_metrics.py to generate BAM summary, gene-level, and per-base coverage metrics.

Command
```bash
python extract_cmv_bam_metrics.py
```
What this script does

This script will generate:
```bash
<sample>_summary.csv
<sample>_gene.csv
<sample>_per_base.csv
```
These files describe:

total reads
average depth
gene-level depth
breadth of coverage
per-base depth
9. Summarise BLAST Hits

Use summarize_blast_hits.py to convert raw BLAST TSV files into simpler CSV outputs.

Command
```bash
python summarize_blast_hits.py --input_dir blast_results/blast_output --output_dir blast_results
```
What this script does

This script will generate:

```bash 
<sample>_blast_top5.csv
<sample>_blast_heatmap.csv
blast_heatmap.csv
```

These outputs make BLAST results easier to review.

10. Launch the Interactive Report

Use app.py to review the pipeline outputs in a Streamlit report.

Command
```bash
streamlit run app.py
``` 
What this report shows

The report can display:

```bash
gene-level coverage
per-base depth
alignment summary metrics
BLAST summaries
raw BLAST hits
variant frequency tables
FastQC reports
annotated VCF results
```

Example Full Workflow

A typical workflow may look like this:

```bash
bash filtering.sh input_fastq/ qc_output/ 20 4
bash alignment_bwa-mem.sh qc_output/filtered_reads reference.fa output/ --threads 4
bash unaligned_reads.sh output/BAM output/ 4
bash 02_blast_cmv_summary.sh input_fastq/ blast_db/ blast_results/ 4
bash variant_calling_vaf.sh output/BAM reference.fa variant_results/ 4
bash annotate_vcfs_snpeff.sh variant_results/vcf snpeff_db/ 4
streamlit run app.py
```

Error Handling

If an error happens during command execution, the script will usually stop and print an error message in the terminal.

Common errors include:

Invalid file paths: check that input folders and files exist
Missing arguments: make sure all required command arguments are given
Missing software: ensure required tools are installed in the environment
Reference problems: confirm the reference FASTA and related databases were prepared correctly
Empty outputs: check that the input FASTQ files contain reads and are in the expected format
Output Files

The pipeline produces several important output files:

*_summary.csv: sample-level summary metrics
*_gene.csv: gene-level depth and breadth
*_per_base.csv: per-base depth values
*_blast.tsv: raw BLAST results
*_blast_top5.csv: top BLAST matches
*_vaf.tsv: variant allele frequency table
*.annot.vcf: annotated variant files
Notes for Users
Run the steps in order unless you already have the required intermediate files
Always check that the reference genome and databases were prepared before starting
Keep output folders organised so results from different runs are not mixed together
Use the Streamlit report at the end to review the outputs more easily