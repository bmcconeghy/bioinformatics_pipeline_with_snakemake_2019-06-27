# Building a bioinformatics pipeline with snakemake
Introduction to snakemake for processing biological sequencing data.
## About the lesson
This course is an introduction to using the Snakemake workflow management system to create reproducible and scalable bioinformatics pipelines. Workflows are written in Python and can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. In this course, we will be writing a Snakemake pipeline that takes DNA reads as input and runs basic QC on them, maps them to a genome of interest, sorts them, indexes them, and calls genomic variants.

Most of this course is based on the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) found on snakemake's readthedocs page.

## Prerequisites
Familiarity with Python, Bash, and some bioinformatics tools (cutadapt, bwa, samtools, etc.).
Software installed in conda environment:
* Python ≥3.6
* Snakemake 5.4.2
* BWA 0.7.12
* SAMtools 1.9
* Pysam 0.15.0
* BCFtools 1.9
* Graphviz 2.38.0
* Jinja2 2.10
* NetworkX 2.1
* Matplotlib 2.2.3

## Getting Started (Cedar)
Installation of miniconda3 is necessary to create the virtual environment we will be working in.
1. Log into the Cedar cluster using the provided credentials.
2. Run the following on the command line: `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`.
3. Once the download has finished, run: `bash Miniconda3-latest-Linux-x86_64.sh` and follow the on-screen instructions. Press `enter` when asked where to install miniconda3 (or specify location).

## Getting Started (Local)
If Cedar is not working, one can follow along by running this locally as well. For Mac OS X, skip this section, for Windows 10 users, do the following:
1. Follow the instructions on [this website](https://itsfoss.com/install-bash-on-windows/).
2. Essentially, what the above website will say is to first enable the linux subsystem by running this command: `Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux` in a administrator PowerShell session.
3. Restart your computer.
4. Install Ubuntu from the Microsoft store.
5. Follow steps 2 and 3 from 'Getting Started (Cedar)'.
6. Say yes to conda initialization and open a new session for the settings to apply.

## Setting up snakemake
1. Clone this workshop's repo: `git clone https://github.com/bmcconeghy/bioinformatics_pipeline_with_snakemake_2019-06-27.git`.
2. Change to newly created directory: `cd bioinformatics_pipeline_with_snakemake_2019-06-27`.
3. Add the bioconda and conda-forge channels: `conda config --add channels bioconda` and `conda --add channels conda-forge`.
4. Use the environment file to download and install all necessary packages for this workshop into a conda environment: `conda create -f smk_542_env.yaml`. This may take a few minutes (total download size is ~265MB).
4. a. If the previous step did not work, run: `conda create -n smk_542 snakemake=5.4.2 bwa=0.7.12 samtools=1.9 pysam=0.15.0 bcftools=1.9 graphviz=2.38.0 jinja2=2.10 networkx=2.1 matplotlib=2.2.3`
5. Once the installation is complete, activate the environment: `conda activate smk_542`. You now have access to all the packages installed!
6. Unzip the tarball: `tar -xzvf data.tar.gz`
7. The data we will be using is stored in the `data` directory.

## Step 1: Mapping reads
1. We will begin by writing a rule to map Illumina reads (single-end) to a reference genome.
2. To do so, we will use the widely used tool, BWA, specifically its MEM subcommand.
3. In the working directory (snakmake-tutorial), create a new file called `Snakefile`. I like using VS Code as an editor, but you may use what you wish. On UNIX servers though, `vim` is ubiquitous.
4. In the Snakefile, define the following rule:
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```
5. The above code is a basic version of a Snakemake rule. It has a name, inputs, ouptput, and a shell directive.
6. When a workflow is executed, Snakemake tries to generate the given target files. These can be specified via the command line. Try running: `snakemake -np mapped_reads/A.bam`
7. This produces the snakemake output as if it were run, but without actually running it; a so called 'dry-run'.
8. Now, re-run the above command without the `-np` parameters: `snakemake mapped_reads/A.bam`
9. 
