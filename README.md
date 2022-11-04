
# kGWASflow    (🚧🚧  Under Development 🚧🚧) <img align="right" width="300" src="https://user-images.githubusercontent.com/42179487/194161153-cc832e57-dd03-481b-8eed-34cb13ba3097.png">

A modular, flexible and reproducible Snakemake workflow to perform k-mers-based GWAS.

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.14-blue.svg)](https://snakemake.github.io) ![License](https://img.shields.io/badge/license-MIT-blue.svg)
[![GitHub actions status](https://github.com/akcorut/kGWASflow/workflows/Tests/badge.svg?branch=main)](https://github.com/akcorut/kGWASflow/actions)
[![DOI](https://zenodo.org/badge/421139649.svg)](https://zenodo.org/badge/latestdoi/421139649)

## Table of Contents

* [Summary](#summary)
* [Installation](#installation)
* [Configuration](#configuration)
* [Usage](#usage)
* [Authors](#authors)
* [License](#license)

## Summary

**kGWASflow** is a [Snakemake](https://snakemake.github.io) pipeline developed for performing k-mers-based genome-wide association study (GWAS) based on the method developed by [Voichek et al. (2020)](https://www.nature.com/articles/s41588-020-0612-7). It performs several pre-GWAS analysis including read trimming, quality control and k-mer counting. It implements the [kmersGWAS method worfklow](https://github.com/voichek/kmersGWAS/blob/master/manual.pdf) for performing k-mers-based GWAS. The pipeline also contains post-GWAS analysis, such as mapping k-mers to a reference genome, finding and mapping the source reads for k-mers, assembling source reads of k-mers into contigs and mapping them to a reference genome. kGWASflow is also highly customizable and offers users multiple options to choose from depends on their needs.

![My project-1-3](https://user-images.githubusercontent.com/42179487/198741149-406abb40-5d1c-4ed0-9a2f-1c6fd9ebed3c.png)

___________

## Installation

### Step 1: Obtain the latest release of this workflow

**1. Clone this repository to your local machine using below command:**

```bash
git clone https://github.com/akcorut/kGWASflow.git
```

**Alternatively**, you can also download and extract the [source code of the latest release](https://github.com/akcorut/kGWASflow/releases).

**2. Change into the kGWASflow directory:**

```bash
cd kGWASflow
```

### Step 2: Install Snakemake and the other dependencies

In order to use this worklow, you need [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to be installed (to install `conda`, please follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)).

**1. Install Snakemake via [`mamba`](https://github.com/mamba-org/mamba) package manager:**

Snakemake recommends `mamba` to be used to install snakemake. More detailed information can be found in the [Snakemake manual](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). To install `mamba`, you can use the below command:

```bash
conda install -c conda-forge mamba
```

After installing `mamba`, you can use below commands to install and activate snakemake and the other dependencies:

```bash
mamba env create -f environment.yaml
conda activate kGWASflow
```

**1a. Alternative installation without `mamba`** 

You can also install snakemake and the other dependencies without `mamba` and just using `conda` as below:

```bash
# This assumes conda is installed in your local machine or computing environment
conda env create -f environment.yaml
conda activate kGWASflow
```

### Other Options: 

The other options on how to deploy this workflow can be found in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=akcorut%2FkGWASflow).

___________

## Configuration

Configure the workflow according to your needs by modifying the files in the `config/` folder.

- `config/config.yaml` is a YAML file containing the workflow configuration.

- `config/samples.tsv` is a TSV file containing the sample information.

- `config/phenos.tsv` is a TSV file contains the phenotype information.

For more information, please click [here](https://github.com/akcorut/kGWASflow/tree/main/config#configuration-settings).

___________

## Usage

After [changing into the kGWASflow directory](https://github.com/akcorut/kGWASflow#step-1-obtain-the-latest-release-of-this-workflow) and [activating the kGWASflow conda environment](https://github.com/akcorut/kGWASflow#step-2-install-snakemake-and-the-other-dependencies), you can start using the workflow as below:

**1. Test your configuration by performing a dry-run**

```bash
snakemake -n --use-conda 
```

**2. Run the workflow and install software dependencies**

```bash
snakemake --cores all --use-conda
```

If you want to run the workflow with a different `config.yaml` file, you can us the `--configfile` parameter to specify it:

```bash
snakemake --use-conda --configfile <path/to/config.yaml>
```

The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=akcorut%2FkGWASflow).

___________

## Authors

kGWASflow was developed by [Adnan Kivanc Corut ](https://www.github.com/akcorut).

___________

## Citation

If you use kGWASflow in your research, please cite using the DOI: [10.5281/zenodo.5790132](https://doi.org/10.5281/zenodo.5790132) and the original method paper by [Voichek et al. (2020)](https://www.nature.com/articles/s41588-020-0612-7):

> Kivanc Corut. (2021). akcorut/kGWASflow: Version 0.1.0-beta (Pre-release) (v0.1.0-beta). Zenodo. https://doi.org/10.5281/zenodo.5790132

> Voichek, Y., Weigel, D. Identifying genetic variants underlying phenotypic variation in plants without complete genomes.  
> Nat Genet 52, 534–540 (2020). https://doi.org/10.1038/s41588-020-0612-7

## License
kGWASflow is licensed under the [MIT](LICENSE.md) license.
