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
* FastQC 0.11.8
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
4. a. If the previous step did not work, run: `conda create -n smk_542 snakemake=5.4.2 bwa=0.7.12 samtools=1.9 pysam=0.15.0 bcftools=1.9 fastqc=0.11.8 graphviz=2.38.0 jinja2=2.10 networkx=2.1 matplotlib=2.2.3`
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
8. One thing you'll notice is some output from BWA, like:
```
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.672 CPU sec, 0.670 real sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/A.fastq
[main] Real time: 0.790 sec; CPU: 0.719 sec
```
You can deal with this neatly by adding a `log` directive to the `bwa_map` rule:
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    log:
        "log/bwa_map/{sample}.log"
    shell:
        "bwa mem {input} 2> {log} | samtools view -Sb - > {output}"
```
9. Now, when you actually run snakemake (without `-n`), the stderr messages are piped into the log file (for each sample requested).

## Step 3: QC of reads
It is always a good idea to take a look at the quality of the reads you will be analyzing. A common tool to do so is FastQC. Normally, this would be the first step of a bioinformatics analysis pipeline that deals with NGS reads, but I wanted to go through some of Snakemake's concepts first, so it's listed here as step 3.
1. Let's put this rule at the top of the Snakefile:
```
rule fastqc:
    input:
        "data/samples/{sample}.fastq"
    output:
        "fastqc/{sample}_fastqc.html",
        "fastqc/{sample}_fastqc.zip"
    params:
        "--outdir fastqc"
    log:
        "log/fastqc/{sample}.log"
    shell:
        "fastqc {params} {input} &> {log}"
```
2. This rule has many of the directives we have seen so far, as well as the `params` one.
3. The `params` directive is essentially a way to make your code more readable and flexible for tweaking down the road.
4. Note how `output` does not even have to be specified in the `shell` directive.
5. Try creating the fastqc output for the `A.fastq` sample. Hint: You only need to specify one of the outputs as both will be created regardless (this is how the tool runs by default).

## Step 4: Sorting read alignments
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
4. Try running snakemake, asking for the sorted version of B.bam: `snakemake -np sorted_reads/B.bam`. Notice how snakemake knows it needs to run `bwa_map` first in order to create `mapped_reads/B.bam`?

## Step 5: Indexing read alignments and visualizing the DAG of jobs
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
![DAG with 3 rules](https://github.com/bmcconeghy/bioinformatics_pipeline_with_snakemake_2019-06-27/blob/master/examples/dag_3rules.svg)
5. Notice that samples for which rules have been run have a dotted line as their outline.

## Step 6: Calling genomic variants
This next step will aggregate the mapped reads from all samples and jointly call genomic variants on them. This will use samtools and bcftools.
1. We will use a built in function of snakemake called `expand`, in conjunction with a list of samples (instead of generalized wilcards). This is similar to the functionality of formatted strings in python.
2. At the very top of your Snakefile, add the list: `SAMPLES = ["A", "B"]`.
3. Now that we have defined what samples will be run through the pipeline, we can write the `bcftools_call` rule:
```
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```
4. Notice that this rule has **named** inputs and how they get used in the `shell` directive. Also, note that you can combine named and not named input and output files in the same rule (although the positional ones must come first; equivalent to python functions with keyword arguments).
5. See if you can produce a DAG that is equivalent to this one:
![DAG with bcftools](https://github.com/bmcconeghy/bioinformatics_pipeline_with_snakemake_2019-06-27/blob/master/examples/dag_bcftools.svg)

## Step 7: Custom scripts
Usually, a workflow not only consists of invoking various tools, but also contains custom code to e.g. calculate summary statistics or create plots. Although snakemake allows you to run python code directly with the `run` directive (instead of `shell`), it is  usually reasonable to move such logic into separate scripts. For this, we will use the `script` directive.
1. Add the following rule to your Snakefile:
```
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```
2. This rule will generate a histogram of the quality scores that have been assigned to the variant calls in the `calls/all.vcf` file.
3. The actual script will be stored in a subdirectory called `scripts`. Create this directory: `mkdir scripts`. Note: script paths are always relative to the Snakefile.
4. Create the file `plot-quals.py` in the `scripts` directory with the following code:
```
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])
```
5. Although there are other strategies to invoke separate scripts from your workflow (e.g., invoking them via shell commands), the benefit of this is obvious: the script logic is separated from the workflow logic (and can be even shared between workflows), but boilerplate code like the parsing of command line arguments is unnecessary.
6. You can also run R scripts this way as well, although the syntax is slightly different.

## Step 8: Adding a target rule
So far, we always executed the workflow by specifying a target file at the command line. Apart from filenames, Snakemake also accepts **rule names as targets** if the referred rule does not have wildcards.
1. It is best practice to have a rule `all` at the top of the workflow which hass all typically desired target files as input files.
2. In this specific case, we add the following to the top of the Snakefile:
```
rule all:
    input:
        "plots/quals.svg"
```
3. And of course, to test it, run: `snakemake -n`. Note how you do not have to specify a target explicitly as this is done within the Snakefile.

## Entire workflow
Your Snakefile (and thus, workflow) should look like this:
```
SAMPLES = ["A", "B"]


rule all:
    input:
        "plots/quals.svg"


rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    log:
        "log/bwa_map/{sample}.log"
    shell:
        "bwa mem {input} 2> {log} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```

# Advanced
## Step 1: Number of threads
1. For some tools, the number of threads can be specified to increase the speed of computation.
2. You can specify a directive, `threads`, that you can use 
