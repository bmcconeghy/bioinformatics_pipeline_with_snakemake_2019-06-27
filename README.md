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

## (Option A) Getting Started - Cedar
Installation of miniconda3 is necessary to create the virtual environment we will be working in.
1. Log into the Cedar cluster using the provided credentials.
2. Run the following on the command line: `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`.
3. Once the download has finished, run: `bash Miniconda3-latest-Linux-x86_64.sh` and follow the on-screen instructions. Press `enter` when asked where to install miniconda3 (or specify location).

## (Option B) Getting Started - Local
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
4. Use the environment file to download and install all necessary packages for this workshop into a conda environment: `conda env create -f smk_542_env.yaml`. This may take a few minutes (total download size is ~265MB).
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
8. Now, re-run the above command without the `-np` arguments: `snakemake mapped_reads/A.bam`
9. The shell command invokes bwa mem with reference genome and reads, and pipes the output into samtools, which creates a compressed BAM file containing the alignments. The output of samtools is then piped into the output file defined by the rule.
10. to generate the target files, snakemake applies the rules given in the Snakefile in a top-down way. For each input file of a job, Snakemake again (i.e. recursively) determines rules that can be applied to generate it. This yields a directed acyclic graph (DAG) of jobs where the edges represent dependencies.

## Step 2: Generalizing the mapping rule
1. Snakemake allows for generalizations for the inputs and outputs by way of **wildcards**. These are denoted within curly braces {}.
2. We will now simply replace the `A` with `{sample}` in the Snakefile. So the resulting code will look like:
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```
3. Once snakemake determines that `{sample}` can be replaced with an appropriate value in the output file, it will propagate that value to all occurrences of `{sample}` in the input files, and thereby determine the necessary input for the resulting job.
4. Try running snakemake with this edited Snakefile with `B.bam` as your requested output: `snakemake -np mapped_reads/B.bam`.
5. You can also specify multiple outputs: `snakemake -np mapped_reads/A.bam mapped_reads/B.bam`
6. Note that snakemake only proposes to create the output file `mapped_reads/B.bam` because we already created `A.bam` previously. This is because the output is newer than the input. You can trigger rerunning the alignment for `A.bam` by 'touching' it: `touch data/samples/A.fastq`.
7. Try running  `snakemake -np mapped_reads/A.bam mapped_reads/B.bam` again. You'll see that it will now propose running both alignments.

## Step 3: Sorting read alignments
For downstream analyses (and in many cases), aligned reads need to be sorted. This can achieved with the **samtools** tool.
1. Add another rule below the `bwa_map` rule, to sort the reads. Call it `samtools_sort`:
```
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```
2. This rule will take any bam file from the `mapped_reads` directory and create a sorted version under `sorted_reads`. One great feature of snakemake is that it creates missing directories!
3. The `-T` parameter allows you to specify a prefix, in this case, the `sorted_reads` directory.
4. Try running snakemake, asking for the sorted version of B.bam: `snakemake -np sorted_reads/B.bam`. Notice how it knows it needs to run `bwa_map` first in order to create `mapped_reads/B.bam`?

## Step 4: Indexing read alignments and visualizing the DAG of jobs
Once again, we use samtools to index the sorted reads (for random access).
1. Write a rule below your sorting rule. Call it `samtools_index`:
```
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```
2. At this point in time, since we have 3 steps in our pipeline, now would be a great time to create a DAG to visualize what is happening in a graphical way: `snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg`.
3. The above command uses Graphviz's dot command to process the ouput from snakemake (using its `--dag` argument) into an svg file.
4. Once you've run this command, you can open this file to take a look. It should look something like this:
![alt text](data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhRE9DVFlQRSBzdmcgUFVCTElDICItLy9XM0MvL0RURCBTVkcgMS4xLy9FTiIKICJodHRwOi8vd3d3LnczLm9yZy9HcmFwaGljcy9TVkcvMS4xL0RURC9zdmcxMS5kdGQiPgo8IS0tIEdlbmVyYXRlZCBieSBncmFwaHZpeiB2ZXJzaW9uIDIuMzguMCAoMjAxNDA0MTMuMjA0MSkKIC0tPgo8IS0tIFRpdGxlOiBzbmFrZW1ha2VfZGFnIFBhZ2VzOiAxIC0tPgo8c3ZnIHdpZHRoPSIyNzBwdCIgaGVpZ2h0PSIxODhwdCIKIHZpZXdCb3g9IjAuMDAgMC4wMCAyNzAuMDAgMTg4LjAwIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIj4KPGcgaWQ9ImdyYXBoMCIgY2xhc3M9ImdyYXBoIiB0cmFuc2Zvcm09InNjYWxlKDEgMSkgcm90YXRlKDApIHRyYW5zbGF0ZSg0IDE4NCkiPgo8dGl0bGU+c25ha2VtYWtlX2RhZzwvdGl0bGU+Cjxwb2x5Z29uIGZpbGw9IndoaXRlIiBzdHJva2U9Im5vbmUiIHBvaW50cz0iLTQsNCAtNCwtMTg0IDI2NiwtMTg0IDI2Niw0IC00LDQiLz4KPCEtLSAwIC0tPgo8ZyBpZD0ibm9kZTEiIGNsYXNzPSJub2RlIj48dGl0bGU+MDwvdGl0bGU+CjxwYXRoIGZpbGw9Im5vbmUiIHN0cm9rZT0iI2FmZDg1NiIgc3Ryb2tlLXdpZHRoPSIyIiBkPSJNMTEwLC0zNkMxMTAsLTM2IDEyLC0zNiAxMiwtMzYgNiwtMzYgLTcuMTA1NDNlLTE1LC0zMCAtNy4xMDU0M2UtMTUsLTI0IC03LjEwNTQzZS0xNSwtMjQgLTcuMTA1NDNlLTE1LC0xMiAtNy4xMDU0M2UtMTUsLTEyIC03LjEwNTQzZS0xNSwtNiA2LC0wIDEyLC0wIDEyLC0wIDExMCwtMCAxMTAsLTAgMTE2LC0wIDEyMiwtNiAxMjIsLTEyIDEyMiwtMTIgMTIyLC0yNCAxMjIsLTI0IDEyMiwtMzAgMTE2LC0zNiAxMTAsLTM2Ii8+Cjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIHg9IjYxIiB5PSItMTUuNSIgZm9udC1mYW1pbHk9InNhbnMiIGZvbnQtc2l6ZT0iMTAuMDAiPnNhbXRvb2xzX2luZGV4PC90ZXh0Pgo8L2c+CjwhLS0gMSAtLT4KPGcgaWQ9Im5vZGUyIiBjbGFzcz0ibm9kZSI+PHRpdGxlPjE8L3RpdGxlPgo8cGF0aCBmaWxsPSJub25lIiBzdHJva2U9IiNhZmQ4NTYiIHN0cm9rZS13aWR0aD0iMiIgZD0iTTI1MCwtMzZDMjUwLC0zNiAxNTIsLTM2IDE1MiwtMzYgMTQ2LC0zNiAxNDAsLTMwIDE0MCwtMjQgMTQwLC0yNCAxNDAsLTEyIDE0MCwtMTIgMTQwLC02IDE0NiwtMCAxNTIsLTAgMTUyLC0wIDI1MCwtMCAyNTAsLTAgMjU2LC0wIDI2MiwtNiAyNjIsLTEyIDI2MiwtMTIgMjYyLC0yNCAyNjIsLTI0IDI2MiwtMzAgMjU2LC0zNiAyNTAsLTM2Ii8+Cjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIHg9IjIwMSIgeT0iLTE1LjUiIGZvbnQtZmFtaWx5PSJzYW5zIiBmb250LXNpemU9IjEwLjAwIj5zYW10b29sc19pbmRleDwvdGV4dD4KPC9nPgo8IS0tIDIgLS0+CjxnIGlkPSJub2RlMyIgY2xhc3M9Im5vZGUiPjx0aXRsZT4yPC90aXRsZT4KPHBhdGggZmlsbD0ibm9uZSIgc3Ryb2tlPSIjZDg1NjU2IiBzdHJva2Utd2lkdGg9IjIiIGQ9Ik0xMDYsLTEwOEMxMDYsLTEwOCAxNiwtMTA4IDE2LC0xMDggMTAsLTEwOCA0LC0xMDIgNCwtOTYgNCwtOTYgNCwtODQgNCwtODQgNCwtNzggMTAsLTcyIDE2LC03MiAxNiwtNzIgMTA2LC03MiAxMDYsLTcyIDExMiwtNzIgMTE4LC03OCAxMTgsLTg0IDExOCwtODQgMTE4LC05NiAxMTgsLTk2IDExOCwtMTAyIDExMiwtMTA4IDEwNiwtMTA4Ii8+Cjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIHg9IjYxIiB5PSItODcuNSIgZm9udC1mYW1pbHk9InNhbnMiIGZvbnQtc2l6ZT0iMTAuMDAiPnNhbXRvb2xzX3NvcnQ8L3RleHQ+CjwvZz4KPCEtLSAyJiM0NTsmZ3Q7MCAtLT4KPGcgaWQ9ImVkZ2UxIiBjbGFzcz0iZWRnZSI+PHRpdGxlPjImIzQ1OyZndDswPC90aXRsZT4KPHBhdGggZmlsbD0ibm9uZSIgc3Ryb2tlPSJncmV5IiBzdHJva2Utd2lkdGg9IjIiIGQ9Ik02MSwtNzEuNjk2NkM2MSwtNjMuOTgyNyA2MSwtNTQuNzEyNSA2MSwtNDYuMTEyNCIvPgo8cG9seWdvbiBmaWxsPSJncmV5IiBzdHJva2U9ImdyZXkiIHN0cm9rZS13aWR0aD0iMiIgcG9pbnRzPSI2NC41MDAxLC00Ni4xMDQzIDYxLC0zNi4xMDQzIDU3LjUwMDEsLTQ2LjEwNDQgNjQuNTAwMSwtNDYuMTA0MyIvPgo8L2c+CjwhLS0gMyAtLT4KPGcgaWQ9Im5vZGU0IiBjbGFzcz0ibm9kZSI+PHRpdGxlPjM8L3RpdGxlPgo8cGF0aCBmaWxsPSJub25lIiBzdHJva2U9IiNkODU2NTYiIHN0cm9rZS13aWR0aD0iMiIgZD0iTTI0NiwtMTA4QzI0NiwtMTA4IDE1NiwtMTA4IDE1NiwtMTA4IDE1MCwtMTA4IDE0NCwtMTAyIDE0NCwtOTYgMTQ0LC05NiAxNDQsLTg0IDE0NCwtODQgMTQ0LC03OCAxNTAsLTcyIDE1NiwtNzIgMTU2LC03MiAyNDYsLTcyIDI0NiwtNzIgMjUyLC03MiAyNTgsLTc4IDI1OCwtODQgMjU4LC04NCAyNTgsLTk2IDI1OCwtOTYgMjU4LC0xMDIgMjUyLC0xMDggMjQ2LC0xMDgiLz4KPHRleHQgdGV4dC1hbmNob3I9Im1pZGRsZSIgeD0iMjAxIiB5PSItODcuNSIgZm9udC1mYW1pbHk9InNhbnMiIGZvbnQtc2l6ZT0iMTAuMDAiPnNhbXRvb2xzX3NvcnQ8L3RleHQ+CjwvZz4KPCEtLSAzJiM0NTsmZ3Q7MSAtLT4KPGcgaWQ9ImVkZ2UyIiBjbGFzcz0iZWRnZSI+PHRpdGxlPjMmIzQ1OyZndDsxPC90aXRsZT4KPHBhdGggZmlsbD0ibm9uZSIgc3Ryb2tlPSJncmV5IiBzdHJva2Utd2lkdGg9IjIiIGQ9Ik0yMDEsLTcxLjY5NjZDMjAxLC02My45ODI3IDIwMSwtNTQuNzEyNSAyMDEsLTQ2LjExMjQiLz4KPHBvbHlnb24gZmlsbD0iZ3JleSIgc3Ryb2tlPSJncmV5IiBzdHJva2Utd2lkdGg9IjIiIHBvaW50cz0iMjA0LjUsLTQ2LjEwNDMgMjAxLC0zNi4xMDQzIDE5Ny41LC00Ni4xMDQ0IDIwNC41LC00Ni4xMDQzIi8+CjwvZz4KPCEtLSA0IC0tPgo8ZyBpZD0ibm9kZTUiIGNsYXNzPSJub2RlIj48dGl0bGU+NDwvdGl0bGU+CjxwYXRoIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzU2ZDhhOSIgc3Ryb2tlLXdpZHRoPSIyIiBkPSJNOTEsLTE4MEM5MSwtMTgwIDMxLC0xODAgMzEsLTE4MCAyNSwtMTgwIDE5LC0xNzQgMTksLTE2OCAxOSwtMTY4IDE5LC0xNTYgMTksLTE1NiAxOSwtMTUwIDI1LC0xNDQgMzEsLTE0NCAzMSwtMTQ0IDkxLC0xNDQgOTEsLTE0NCA5NywtMTQ0IDEwMywtMTUwIDEwMywtMTU2IDEwMywtMTU2IDEwMywtMTY4IDEwMywtMTY4IDEwMywtMTc0IDk3LC0xODAgOTEsLTE4MCIvPgo8dGV4dCB0ZXh0LWFuY2hvcj0ibWlkZGxlIiB4PSI2MSIgeT0iLTE2NSIgZm9udC1mYW1pbHk9InNhbnMiIGZvbnQtc2l6ZT0iMTAuMDAiPmJ3YV9tYXA8L3RleHQ+Cjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIHg9IjYxIiB5PSItMTU0IiBmb250LWZhbWlseT0ic2FucyIgZm9udC1zaXplPSIxMC4wMCI+c2FtcGxlOiBBPC90ZXh0Pgo8L2c+CjwhLS0gNCYjNDU7Jmd0OzIgLS0+CjxnIGlkPSJlZGdlMyIgY2xhc3M9ImVkZ2UiPjx0aXRsZT40JiM0NTsmZ3Q7MjwvdGl0bGU+CjxwYXRoIGZpbGw9Im5vbmUiIHN0cm9rZT0iZ3JleSIgc3Ryb2tlLXdpZHRoPSIyIiBkPSJNNjEsLTE0My42OTdDNjEsLTEzNS45ODMgNjEsLTEyNi43MTIgNjEsLTExOC4xMTIiLz4KPHBvbHlnb24gZmlsbD0iZ3JleSIgc3Ryb2tlPSJncmV5IiBzdHJva2Utd2lkdGg9IjIiIHBvaW50cz0iNjQuNTAwMSwtMTE4LjEwNCA2MSwtMTA4LjEwNCA1Ny41MDAxLC0xMTguMTA0IDY0LjUwMDEsLTExOC4xMDQiLz4KPC9nPgo8IS0tIDUgLS0+CjxnIGlkPSJub2RlNiIgY2xhc3M9Im5vZGUiPjx0aXRsZT41PC90aXRsZT4KPHBhdGggZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTZkOGE5IiBzdHJva2Utd2lkdGg9IjIiIGQ9Ik0yMzEsLTE4MEMyMzEsLTE4MCAxNzEsLTE4MCAxNzEsLTE4MCAxNjUsLTE4MCAxNTksLTE3NCAxNTksLTE2OCAxNTksLTE2OCAxNTksLTE1NiAxNTksLTE1NiAxNTksLTE1MCAxNjUsLTE0NCAxNzEsLTE0NCAxNzEsLTE0NCAyMzEsLTE0NCAyMzEsLTE0NCAyMzcsLTE0NCAyNDMsLTE1MCAyNDMsLTE1NiAyNDMsLTE1NiAyNDMsLTE2OCAyNDMsLTE2OCAyNDMsLTE3NCAyMzcsLTE4MCAyMzEsLTE4MCIvPgo8dGV4dCB0ZXh0LWFuY2hvcj0ibWlkZGxlIiB4PSIyMDEiIHk9Ii0xNjUiIGZvbnQtZmFtaWx5PSJzYW5zIiBmb250LXNpemU9IjEwLjAwIj5id2FfbWFwPC90ZXh0Pgo8dGV4dCB0ZXh0LWFuY2hvcj0ibWlkZGxlIiB4PSIyMDEiIHk9Ii0xNTQiIGZvbnQtZmFtaWx5PSJzYW5zIiBmb250LXNpemU9IjEwLjAwIj5zYW1wbGU6IEI8L3RleHQ+CjwvZz4KPCEtLSA1JiM0NTsmZ3Q7MyAtLT4KPGcgaWQ9ImVkZ2U0IiBjbGFzcz0iZWRnZSI+PHRpdGxlPjUmIzQ1OyZndDszPC90aXRsZT4KPHBhdGggZmlsbD0ibm9uZSIgc3Ryb2tlPSJncmV5IiBzdHJva2Utd2lkdGg9IjIiIGQ9Ik0yMDEsLTE0My42OTdDMjAxLC0xMzUuOTgzIDIwMSwtMTI2LjcxMiAyMDEsLTExOC4xMTIiLz4KPHBvbHlnb24gZmlsbD0iZ3JleSIgc3Ryb2tlPSJncmV5IiBzdHJva2Utd2lkdGg9IjIiIHBvaW50cz0iMjA0LjUsLTExOC4xMDQgMjAxLC0xMDguMTA0IDE5Ny41LC0xMTguMTA0IDIwNC41LC0xMTguMTA0Ii8+CjwvZz4KPC9nPgo8L3N2Zz4K "DAG with 3 rules")
