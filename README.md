# Building a bioinformatics pipeline with snakemake
Introduction to snakemake for processing biological sequencing data.
## About the lesson
This course is an introduction to using the Snakemake workflow management system to create reproducible and scalable bioinformatics pipelines. Workflows are written in Python and can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. In this course, we will be writing a Snakemake pipeline that takes DNA reads as input and runs basic QC on them, maps them to a genome of interest, sorts them, indexes them, and calls genomic variants.

Most of this course is based on the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) found on snakemake's readthedocs page

## Prerequisites
Familiarity with Python, Bash, and some bioinformatics tools (cutadapt, bwa, samtools, etc.).
Software:
*Python ≥3.5
*Snakemake 5.2.3
*BWA 0.7.12
*SAMtools 1.3.1
*Pysam 0.15.0
*BCFtools 1.3.1
*Graphviz 2.38.0
*Jinja2 2.10
*NetworkX 2.1
*Matplotlib 2.2.3

## Getting Started
Installation of miniconda3 is necessary to create the virtual environment we will be working in.
1. Log into the Cedar (or Graham) cluster using the provided credentials.
2. Run the following on the command line: `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
3. Once the download has finished, run: `bash Miniconda3-latest-Linux-x86_64.sh` and follow the on-screen instructions. Press `enter` when asked where to install miniconda3 (or specify location).

## Setting up snakemake
4. Create a new directory under your home directory: `mkdir snakemake-tutorial`
5. Change to newly created directory: `cd snakemake-tutorial`
6. 
