# Snakemake workflow: Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This repository contains snakemake workflows for analysing epigenetic signals captured by fragmentation patterns of cell-free DNA.

The following analyses are based on the following publication:

Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises
an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell. 2016 Jan
14;164(1-2):57-68. doi: 10.1016/j.cell.2015.11.050. PubMed PMID: [26771485](http://www.ncbi.nlm.nih.gov/pubmed/26771485)

More details on analyses in this [repository](https://github.com/shendurelab/cfDNA).

*All scripts and binaries are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are __not__ providing any support for these scripts.*

## Table of Contents

- [Snakemake workflow: Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA](#snakemake-workflow-analysis-of-epigenetic-signals-captured-by-fragmentation-patterns-of-cell-free-dna)
    - [Table of Contents](#table-of-contents)
    - [Authors](#authors)
    - [Versions and dependencies](#versions-and-dependencies)
        - [Samtools version used](#samtools-version-used)
    - [Usage](#usage)
        - [Step 1: Obtain a copy of this workflow](#step-1-obtain-a-copy-of-this-workflow)
        - [Step 2: Configure workflow](#step-2-configure-workflow)
        - [Step 3: Install Snakemake](#step-3-install-snakemake)
        - [Step 4: Execute workflow](#step-4-execute-workflow)
    - [Workflows](#workflows)
        - [Windowed protection score overlays](#windowed-protection-score-overlays)
            - [Description](#description)
            - [Input](#input)
            - [Output](#output)
        - [Gene Expression Analysis](#gene-expression-analysis)
            - [Description](#description-1)
            - [Input](#input-1)
            - [Output](#output-1)
        - [Unsupervised analysis](#unsupervised-analysis)
            - [Description](#description-2)
            - [Input](#input-2)
            - [Output](#output-2)

## Authors

- Sebastian Röner (@sroener)

## Versions and dependencies

Our scripts largely depend on Python 2, the [pysam](https://pysam.readthedocs.io/en/latest/api.html), [bx](https://pypi.org/project/bx-python/) and [numpy](https://numpy.org/) libraries. Here a list of versions that we used:

```bash
dependencies:
  - python=2.7.15
  - numpy=1.7.2
  - pysam=0.7.6
  - bx-python=0.8.9
```

We were also using [R 4.0.2](https://cran.r-project.org/), [UCSC binaries](http://hgdownload.cse.ucsc.edu/admin/exe/) for working with bigWig and bigBed files and tools like tabix and samtools from [HTSlib](http://www.htslib.org/).

### Samtools version used

Please note that a samtools binary is included with these scripts. Among other things, this samtools binary allows filtering reads based on insert size/read length. This is an early version of the samtools branch released on [https://github.com/mpieva/samtools-patched](https://github.com/mpieva/samtools-patched). For samtool calls that are not filtering for read length/insert size, other (more recent) versions of samtools might be used. Please note though that parameters might be named differently and that the ASCII encoding of read filters is also special to the samtools version that we used.

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) the repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, `samples.tsv` to specify your sample setup and `regions.tsv` to specify target regions.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

The workflows are executed from the repository root folder. The different analyses have to be executed separately. To specify the respective workflow use the `-s` switch followed by the path of the `Snakefile` (e.g.: `./snakefile_WPS.smk`)

Activate the conda environment:

```bash
conda activate snakemake
```

Test your configuration by performing a dry-run via

```bash
snakemake -s path/to/Snakefile --use-conda -n
```

Execute the workflow locally via

```bash
snakemake -s path/to/Snakefile --use-conda --cores $N
```

using `$N` cores or run it in a cluster environment via

```bash
snakemake -s path/to/Snakefile --use-conda --cluster qsub --jobs 100
```

or

```bash
snakemake -s path/to/Snakefile --use-conda --drmaa --jobs 100
```

If you not only want to fix the software stack but also the underlying OS, use

```bash
snakemake --use-conda --use-singularity
```

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for further details.

## Workflows

### Windowed protection score overlays

#### Description

Fragment endpoint coordinates were extracted from BAM files with the pysam API.Both outer alignment coordinates of PE data were extracted for properly paired reads. Both end coordinates of SR alignments were extracted when PE data was collapsed to SR data by adapter trimming. A fragment’s coverage is defined as all positions between the two (inferred) fragment ends, inclusive of endpoints. We define the Windowed Protection Score (WPS) of a window of size k as the number of molecules spanning the window minus those with an endpoint within the window. We assign the determined WPS to the center of the window.

Windowed protection scores are calculated for all provided regions with additional 1000 bp up- and downstream. The values are then averaged over all provided regions and normalized by substracting the trimmed mean (excluding the outer 25% on both sides of the distribution).

#### Input

- configured by the user ([samples.tsv](config/samples.tsv)):
    - analysis ID
    - samples
    - path to sample .bam files
    - reference samples fro plotting
    - genome build per sample
- configured by the user ([regions.tsv](config/regions.tsv)):
    - bed file containing regions of interest (e.g. TFBS), all having the same length

**Note:** More information about config files [here](config/README.md)

#### Output

- table containing bp specific WPS for regions listed in bed
- line plot showing normalized WPS of multiple samples

### Gene Expression Analysis

The analysis is contained in the `snakefile_GE_analysis.smk` workflow.

#### Description

Human Protein Atlas

FPKM gene expression (GE) values measured for 20,344 Ensembl gene identifiers in 44 human cell lines and 32 primary tissues by the Human Protein Atlas (Uhlén et al., 2015) was downloaded from [http://www.proteinatlas.org/download/rna.csv.zip](http://www.proteinatlas.org/download/rna.csv.zip). Genes with 3 or more non-zero expression values were retained (n=19,378 genes). The GE data set is provided with one decimal precision for the FPKM values. Thus, a zero GE value (0.0) indicates expression in the interval [0, 0.05) Unless otherwise noted, we set the minimum GE value to 0.04 FPKM before log2-transformation.

Fast Fourier transformation (FFT) and smoothing of trajectories

We use parameters to smooth (3 bp Daniell smoother; moving average giving half weight to the end values) and de-trend the data (i.e. subtract the mean of the series and remove a linear trend). A recursive time series filter implemented in R was used to remove high frequency variation from trajectories. 24 filter frequencies (1/seq(5,100,4)) were used, and the first 24 values of the trajectory were taken as init values. The 24-value shift in the resulting trajectories was corrected by repeating the last 24 values of the trajectory.

FFT intensity correlation with expression

WPS was used to calculate periodograms of genomic regions using Fast Fourier Transform (FFT, spec.pgram in R) with frequencies between 1/500 and 1/100 bases. Intensity values for the 120-280 bp frequency range were determined from smooth FFT periodograms. S-shaped Pearson correlation between GE values and FFT intensities was observed around the major inter-nucleosome distance peak, along with a pronounced negative correlation in the 193-199 bp frequency range. The mean intensity in this frequency range was correlated with the average intensity with log2-transformed GE values for downstream analysis.

#### Input

- included in the repository:
    - annotations
    - labels
    - RNAtable from Protein Atlas 
        - blood atlas ["Blood"]
        - protein atlas tissues ["Tissue"]
        - protein atlas tissues + cell lines ["Extended"]
- configured by the user ([samples.tsv](config/samples.tsv)):
    - analysis ID
    - samples
    - path to sample bam files
    - genome build per sample

**Note:** More information about config files [here](config/README.md)

#### Output

- fft_summary tables (results/intermediate/body/fft_summaries)
- plots showing intensities across tissues (results/plots)
- table showing correlation with tissues/cell lines
- table showing correlation rank difference to reference sample

### Unsupervised analysis

#### Description

Workflow containing utility functions to calculate and visualize similarities between samples. 

#### Input

- FFT tables by gene expression analysis

**Note:** FFT tables will be generated automatically if gene expression workflow was not exicuted before.

#### Output

- heatmaps showing spearman correlation between samples from the same experiment (ID) based on configured FFT thresholds
- clustermaps showing spearman correlation between samples from the same experiment (ID) based on configured FFT thresholds with added hierarchical clustering
- Plots showing 2D UMAP projection labeled by kmeans clustering with specified number of clusters
- Plots showing 2D UMAP projection labeled by HDBSCAN clustering