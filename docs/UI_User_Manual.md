# User Interface (UI) Manual for the CMV Whole Genome Sequencing Pipeline

## Table of Contents

- Introduction
- Purpose
- System Requirements
- UI Components Overview
- Common User Interactions
- Troubleshooting UI Issues
- Best Practices
- Feedback and Support

## Introduction

The CMV Whole Genome Sequencing Pipeline user interface provides a desktop-based way to run the pipeline without needing to use the command line for every step. The interface is designed to help users select input files, choose references and databases, run the required analysis steps, monitor progress, and review outputs more easily.

The project also includes an interactive Streamlit report, which allows users to explore the results after the pipeline has finished running.

## Purpose

This manual provides guidance on:

- Navigating the graphical user interface effectively
- Selecting the correct input files and output folders
- Running common pipeline steps from the interface
- Viewing logs and understanding what the pipeline is doing
- Launching the interactive report
- Troubleshooting common user interface issues

## System Requirements

To use the UI effectively, ensure your system meets the following requirements:

- Operating System: Linux recommended
- Python: Python 3.10
- Conda environment: `CMV_Project`
- Required packages installed from `environment.yml`
- Tkinter available for the desktop GUI
- Streamlit installed for the interactive report
- Minimum screen resolution: 1280 x 720 recommended
- Sufficient disk space for FASTQ, BAM, BLAST, and VCF outputs

You will also need the main pipeline dependencies installed, including:

- bwa
- samtools
- fastqc
- trim_galore
- blastn
- seqtk
- lofreq
- snpEff

## UI Components Overview

### 1. Input Selection Section

**Purpose:**  
Allows users to choose the main files and folders needed to run the pipeline.

**Components may include:**

- FASTQ input folder or archive selection
- Reference FASTA selection
- BLAST database selection
- Annotation file selection where needed
- Output directory selection
- Browse buttons for easier path selection

This section is used to tell the software where the input files are and where outputs should be saved.

### 2. Pipeline Steps Section

**Purpose:**  
Allows users to choose which stages of the pipeline to run.

**Examples of pipeline steps available through the GUI include:**

- Quality control and filtering
- Alignment
- Unmapped read extraction
- BLAST analysis
- Variant calling
- Variant annotation
- Summary metric generation

This section helps users run only the steps they need instead of always running the full pipeline.

### 3. Parameters Section

**Purpose:**  
Allows users to enter important settings used during analysis.

**Examples of settings may include:**

- Number of threads
- Quality threshold
- Output folder name
- Port number for the Streamlit report
- Optional compression settings

This section controls how the pipeline runs.

### 4. Log Output Section

**Purpose:**  
Shows live updates while the pipeline is running.

**What it displays:**

- Commands being run
- Progress messages
- Warnings
- Errors
- Completion messages

This is one of the most important areas of the GUI because it helps users see whether the run is progressing normally.

### 5. Action Buttons

**Purpose:**  
Allows users to start tasks and open reports or related tools.

**Examples may include:**

- Run pipeline
- Stop pipeline
- Launch Streamlit report
- Open FastQC report
- Open IGV
- Clear logs

These buttons are used for the main user actions.

### 6. Interactive Report

**Purpose:**  
Provides a visual way to review outputs after the pipeline has run.

**The report can show:**

- Gene-level coverage summaries
- Per-base coverage tables
- Alignment metrics
- BLAST summaries
- Variant allele frequency tables
- Annotated VCF results
- FastQC outputs

This part is usually launched after analysis is complete.

## Common User Interactions

### Selecting Input Data

1. Open the GUI.
2. Locate the input selection section.
3. Choose the folder or archive containing the FASTQ files.
4. Select the CMV reference FASTA.
5. Select the BLAST database folder.
6. Choose the output directory.

Always make sure the selected files match the run you want to perform.

### Running Quality Control and Filtering

1. Select the FASTQ input.
2. Set the quality threshold.
3. Set the number of CPU cores.
4. Select the filtering or QC step.
5. Click the Run button.

The GUI will then run the trimming and FastQC step and display progress in the log window.

### Running Alignment

1. Make sure the reference FASTA has already been indexed.
2. Select the filtered FASTQ folder.
3. Select the alignment step.
4. Set the number of threads.
5. Click Run.

The output BAM files and BAM statistics will be created in the output folder.

### Running BLAST Analysis

1. Select the input FASTQ folder or archive.
2. Select the BLAST database directory.
3. Choose the BLAST step.
4. Set the number of threads.
5. Click Run.

The GUI will generate FASTA files, run BLAST, and create summary CSV outputs.

### Running Variant Calling

1. Select the BAM output folder.
2. Select the CMV reference FASTA.
3. Choose the variant calling step.
4. Set the number of threads.
5. Click Run.

VCF files and VAF tables will be written to the selected output directory.

### Running Variant Annotation

1. Select the VCF directory.
2. Select the snpEff working directory.
3. Choose the annotation step.
4. Click Run.

Annotated VCF files will be produced.

### Launching the Interactive Report

1. Make sure the pipeline outputs have already been generated.
2. Click the button to launch the Streamlit report.
3. If needed, set the port number first.
4. Open the report in the browser.

The report will allow you to inspect the outputs in a more visual way.

### Opening FastQC Results

1. Complete the filtering or FastQC step.
2. Use the button in the GUI to open the FastQC output.
3. Review the HTML report in your browser.

This is useful for checking input read quality.

### Opening IGV

1. Ensure BAM and related output files are available.
2. Click the IGV launch button.
3. Load the relevant files for visual review.

This can be used for manual inspection of alignments and variant evidence.

## Troubleshooting UI Issues

### GUI Does Not Open

**Steps to resolve:**

- Confirm that the conda environment is activated
- Check that Python 3.10 is installed
- Ensure tkinter is available
- Run the script again from the terminal and check for error messages

### Browse Button Not Working

**Steps to resolve:**

- Ensure the system supports file dialog windows
- Check that the selected starting directory exists
- Restart the GUI and try again

### Pipeline Does Not Start

**Steps to resolve:**

- Check that all required input paths have been selected
- Ensure the output directory is writable
- Make sure required software tools are installed in the environment
- Review the log box for error messages

### Logs Show “Tool Not Found”

**Steps to resolve:**

- Confirm that the conda environment is activated
- Reinstall missing dependencies using `environment.yml`
- Test the missing command in the terminal directly, for example:
  - `bwa`
  - `samtools`
  - `blastn`

### Streamlit Report Does Not Open

**Steps to resolve:**

- Make sure `streamlit` is installed
- Check that the selected port is not already in use
- Confirm that output files were successfully created before launching the report
- Try launching manually with:

```bash
streamlit run app.py
```

Results Not Appearing in the Report

Steps to resolve:

Confirm that the output files exist in the selected folders
Check that files such as *_summary.csv, *_gene.csv, *_blast_top5.csv, and *_vaf.tsv were generated
Make sure the correct folder was chosen in the report interface
Refresh or rerun the report after new files are created
FastQC or IGV Buttons Do Not Open Anything

Steps to resolve:

Confirm the related output files exist
Check that the required application is installed
Review the terminal or GUI log for any launch errors
Best Practices
Use clear and organised folder names for each run
Keep input, intermediate, and final output folders separate
Check all selected paths before starting a run
Run the pipeline in the correct order where possible
Review the log window during each run so errors can be seen early
Launch the Streamlit report only after output files are available
Avoid mixing results from multiple runs in the same output directory
Keep your conda environment updated and consistent across runs
Feedback and Support

If you encounter issues or have suggestions for improving the interface:

Review the log output first, as many problems can be identified there
Check that the required software and reference files are installed correctly
Ask the project maintainer or local pipeline administrator for help if needed

When reporting a problem, it is helpful to include:

The step you were trying to run
The exact error message shown in the log
The input files and settings used
Whether the issue happened in the GUI, the report, or both