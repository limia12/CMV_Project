# Installation

These instructions will allow you to install the CMV Whole Genome Sequencing Pipeline and its required dependencies on Linux. The pipeline was developed and tested in a Linux environment, so this is the recommended system for installation and use.

## Pre-requisites

Installation requires:

- Conda
- Python 3.10
- Git
- The project `environment.yml` file
- The project `requirements.txt` file

The conda environment installs the main Python and bioinformatics dependencies used by the pipeline, including:

- bwa
- samtools
- fastqc
- trim_galore
- blastn
- seqtk
- lofreq
- snpEff
- streamlit
- pysam

## Installing

Download the git repository:

```bash
git clone <your_repository_url>
cd CMV_Pipeline
```

Create and activate the conda environment:

```bash 
conda env create -f environment.yml
conda activate CMV_Project
```

Install any additional Python dependencies listed in requirements.txt:

```bash 
pip install -r requirements.txt
``` 

## Reference and Database Setup

Before running the pipeline, you will need to prepare the reference genome and supporting databases.

1. Prepare the reference genome

Run:

```bash
bash index_genome.sh reference.fa
```

This creates the required reference index files for alignment and downstream analysis.

2. Build the snpEff database

Run:

```bash 
bash build_snpeff_db.sh reference.fa annotation.gff snpeff_db/
```

This prepares the local snpEff database used for variant annotation.

3. Create the BLAST database

Run:

```bash
makeblastdb -in viral_sequences.fna -dbtype nucl -out viral_nt_db
```

This creates the local BLAST database used for sequence confirmation.

Final Check

To confirm that installation has worked correctly, check that the main tools are available:

```bash
bwa
samtools
blastn
lofreq
snpEff
streamlit
```

You will also need to ensure that your reference FASTA, annotation GFF, and BLAST database files are available before running the pipeline.