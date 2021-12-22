# kGWASflow    (🚧🚧  Under Development 🚧🚧)

A Snakemake workflow to perform k-mers-based GWAS.

[![DOI](https://zenodo.org/badge/421139649.svg)](https://zenodo.org/badge/latestdoi/421139649)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.10-blue.svg)](https://snakemake.github.io) ![License](https://img.shields.io/badge/license-MIT-blue.svg)
[![Actions Status](https://github.com/akcorut/kGWASflow/workflows/CI/badge.svg)](https://github.com/akcorut/kGWASflow/actions)
[![Actions Status](https://github.com/akcorut/kGWASflow/workflows/Linting/badge.svg)](https://github.com/akcorut/kGWASflow/actions)

## Table of Contents

* [Summary](#Summary)
* [Installation](#installation)
* [Configuration](#configuration)
* [Usage](#usage)
* [Authors](#authors)
* [Acknowledgements](#acknowledgements)
* [License](#license)

## Summary

**kGWASflow** is a [Snakemake](https://snakemake.github.io) pipeline developed for performing k-mers-based genome-wide association study (GWAS) based on the method developed by [Voichek et al. (2020)](https://www.nature.com/articles/s41588-020-0612-7). It performs several pre-GWAS analysis including read trimming, quality control and k-mer counting. It implements the [kmersGWAS method worfklow](https://github.com/voichek/kmersGWAS/blob/master/manual.pdf) for performing k-mers-based GWAS. The pipeline also contains post-GWAS analysis, such as finding and aligning the source reads for k-mers, aligning kmers and source reads to a reference genome.

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

You can also install snakemake and the other dependencies without mamba as below:

```bash
conda env create -f environment.yaml
conda activate kGWASflow
```

## Configuration

## Usage

## Authors

## Acknowledgements

## License
kGWASflow is licensed under the [MIT](LICENSE.md) license.
