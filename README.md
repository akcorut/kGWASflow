
# kGWASflow <img align="right" width="420" src="https://user-images.githubusercontent.com/42179487/194161153-cc832e57-dd03-481b-8eed-34cb13ba3097.png">

A modular, flexible, and reproducible Snakemake workflow to perform k-mers-based GWAS.

[![Anaconda-Server Badge](https://anaconda.org/bioconda/kgwasflow/badges/version.svg)](https://anaconda.org/bioconda/kgwasflow)
![Anaconda-Server Badge](https://anaconda.org/bioconda/kgwasflow/badges/latest_release_date.svg)
![](https://img.shields.io/github/v/tag/akcorut/kGWASflow?label=version&style=flat-square)
[![install with bioconda](https://img.shields.io/badge/Install%20with-conda-brightgreen.svg?style=flat-square)](https://anaconda.org/bioconda/kgwasflow)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/kgwasflow/badges/downloads.svg)](https://anaconda.org/bioconda/kgwasflow)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25-blue.svg)](https://snakemake.github.io)
![License](https://img.shields.io/badge/license-MIT-black.svg)
[![DOI](https://zenodo.org/badge/421139649.svg)](https://zenodo.org/badge/latestdoi/421139649)

## Table of Contents

* [Summary](#summary)
* [Installation](#installation)
* [Configuration](#configuration)
* [Usage](#usage)
* [Testing](#testing)
* [Authors](#authors)
* [Issues](#issues)
* [Maintainers](#maintainers)
* [License](#license)

## Summary

**kGWASflow** is a [Snakemake](https://snakemake.github.io) pipeline developed for performing k-mers-based genome-wide association study (GWAS) based on the method developed by [Voichek et al. (2020)](https://www.nature.com/articles/s41588-020-0612-7). It performs several pre-GWAS analyses, including read trimming, quality control, and k-mer counting. It implements the [kmersGWAS method]([https://github.com/voichek/kmersGWAS/blob/master/manual.pdf](https://www.nature.com/articles/s41588-020-0612-7)) into an easy to use and accessible workflow. The pipeline also contains post-GWAS analyses, such as mapping k-mers to a reference genome, finding and mapping the source reads of k-mers, assembling source reads into contigs, and mapping them to a reference genome. kGWASflow is also highly customizable and offers users multiple options to choose from depending on their needs.

**More information and explanations on how to install, configure and run kGWASflow are provided in the [kGWASflow Wiki](https://github.com/akcorut/kGWASflow/wiki).**

**kGWASflow preprint is out on bioRxiv:** 
- Adnan Kivanc Corut, Jason G. Wallace. kGWASflow: a modular, flexible, and reproducible Snakemake workflow for k-mers-based GWAS (2023). bioRxiv. [https://doi.org/10.1101/2023.07.10.548365](https://doi.org/10.1101/2023.07.10.548365)

![My project-1-3](https://user-images.githubusercontent.com/42179487/198741149-406abb40-5d1c-4ed0-9a2f-1c6fd9ebed3c.png)

___________

## Installation

### Installing via Bioconda 

‼️This is the **preferred method**.‼️

In order to use this workflow, you need [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to be installed (to install `conda`, please follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)).

```bash
# Create a new conda environment with kgwasflow 
# and its dependencies
conda create -c bioconda -n kgwasflow kgwasflow

# Activate kGWASflow conda environment
conda activate kgwasflow

# test kGWASflow
kgwasflow --help
```

### Installing via GitHub (Alternative Method)

Alternatively, kGWASflow can be installed by cloning [the GitHub repository](https://github.com/akcorut/kGWASflow).

```bash
# Clone this repository to your local machin
git clone https://github.com/akcorut/kGWASflow.git

# Change into the kGWASflow directory
cd kGWASflow
```

After cloning the GitHub repo, you can install snakemake and the other dependencies using `conda` as below:

```bash
# This assumes conda is installed
conda env create -f environment.yaml

# Activate kGWASflow conda environment
conda activate kGWASflow
```

Finally, you can install kGWASflow using the setup script as below:

```bash
# Install kgwasflow
python setup.py install

# Test kgwasflow
kgwasflow --help
```

### Other Options: 

The other options on how to deploy this workflow can be found in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=akcorut%2FkGWASflow).

___________

## Configuration

### Initializing kGWASflow

To configure kGWASflow, you first need to initialize a new kGWASflow working directory by following the below steps:

```bash
# Activating the conda environment
conda activate kgwasflow

# Initializing a new kgwasflow working dir
kgwasflow init --work-dir path/to/your/work_dir 
```

or

```bash
# Activating the conda environment
conda activate kgwasflow

# Change into your preferred working directory
cd path/to/your/work_dir 

# Initializing a new kgwasflow working dir
kgwasflow init
```

### kGWASflow Working Directory Structure 

This command will initialize a new kGWASflow working directory with the default configuration files. Below is the directory structure of the working directory:

```
path/to/your/work_dir 
├── config
│   ├── config.yaml
│   ├── phenos.tsv
│   └── samples.tsv
└── test
```

### kGWASflow Configuration Files

Below are the configuration files generated by `kgwasflow init` command:

- `config/config.yaml` is a YAML file containing the workflow configuration.

- `config/samples.tsv` is a TSV file containing the sample information.

- `config/phenos.tsv` is a TSV file containing the phenotype information.

**For more information about each configuration file, please see [kGWASflow Wiki](https://github.com/akcorut/kGWASflow/wiki).**

___________

## Usage

After initializing (`kgwasflow init`) [step](#configuration) and modifying the configuration files, kGWASflow can be run as below:

```shell
# Activating the conda environment
conda activate kgwasflow

# Change into your preferred working directory
cd path/to/your/work_dir 

# Run kgwasflow
kgwasflow run -t 16
```

Below are some of the run examples of kGWASflow:

```shell
Run examples:

1. Run kGWASflow with the default config file, default arguments and 16 threads:

  kgwasflow run -t 16 --snake-default

2. Run kGWASflow with a custom config file and default settings:

  kgwasflow run -t 16 -c path/to/custom_config.yaml

3. Run kGWASflow with user defined working directory:

  kgwasflow run -t 16 --work-dir path/to/work_dir

4. Run kGWASflow in dryrun mode to see what tasks would be executed:

  kgwasflow run -t 16 -n

5. Run kGWASflow using mamba as the conda frontend:

  kgwasflow run -t 16 --conda-frontend mamba

6. Run kGWASflow and generate an HTML report:

  kgwasflow run -t 16 --generate-report
```

Information about how to use kGWASflow with Snakemake commands can be found in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=akcorut%2FkGWASflow).

___________

## Testing

In order to test kGWASflow using an E.coli dataset ([Earle et al. 2016](https://www.nature.com/articles/nmicrobiol201641), [Rahman et al. 2018](https://elifesciences.org/articles/32920))


```shell
# Activating the conda environment
conda activate kgwasflow
```

After activating kGWASflow, you can perform a test run as axplained below:

```shell
Test examples:

1. Run the kGWASflow test in dryrun mode to see what tasks would be executed:

   kgwasflow test -t 16 -n

2. Run the kGWASflow test using the test config file with 16 threads:

   kgwasflow test -t 16

3. Run the kGWASflow test and define the test working directory:

   kgwasflow test -t 16 --work-dir path/to/test_work_dir
```

___________

## Authors

kGWASflow was developed by [Adnan Kivanc Corut](https://www.github.com/akcorut).

___________

## Issues

For Issues: [https://github.com/akcorut/kGWASflow/issues](https://github.com/akcorut/kGWASflow/issues)

**Contributions to the development of kGWASflow are welcome! Create Pull Requests to fix bugs or recommend new features!**

___________

## Maintainers

- [Adnan Kivanc Corut ![orcid]](https://orcid.org/0009-0002-7622-7300), [![github]](https://github.com/akcorut)

[github]: https://cbg-ethz.github.io/V-pipe/img/mark-github.svg
[orcid]: https://cbg-ethz.github.io/V-pipe/img/ORCIDiD_iconvector.svg
___________

## Citation

‼️**kGWASflow preprint is out now on bioRxiv:** [https://doi.org/10.1101/2023.07.10.548365](https://doi.org/10.1101/2023.07.10.548365) ‼️

If you use kGWASflow in your research, please cite using the DOI: [https://doi.org/10.1101/2023.07.10.548365](https://doi.org/10.1101/2023.07.10.548365) and the original method paper by [Voichek et al. (2020)](https://www.nature.com/articles/s41588-020-0612-7):

* Adnan Kivanc Corut, Jason G. Wallace. kGWASflow: a modular, flexible, and reproducible Snakemake workflow for k-mers-based GWAS (2023). bioRxiv. [https://doi.org/10.1101/2023.07.10.548365](https://doi.org/10.1101/2023.07.10.548365)

* Voichek, Y., Weigel, D. Identifying genetic variants underlying phenotypic variation in plants without complete genomes. Nat Genet 52, 534–540 (2020). https://doi.org/10.1038/s41588-020-0612-7

___________

## License
kGWASflow is licensed under the [MIT](LICENSE.md) license.
