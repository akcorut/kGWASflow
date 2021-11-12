# kGWASflow (🚧 Under Development)

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.10-blue.svg)](https://snakemake.github.io) ![License](https://img.shields.io/badge/license-MIT-blue.svg)

## Summary

**kGWASflow** is a [Snakemake](https://snakemake.github.io) pipeline developed for performing k-mers-based genome-wide association study (GWAS) based on the method developed by [Voichek et al. (2020)](https://www.nature.com/articles/s41588-020-0612-7). It performs several pre-GWAS analysis including read trimming, quality control and k-mer counting. It implements the [kmersGWAS method worfklow](https://github.com/voichek/kmersGWAS/blob/master/manual.pdf) for performing k-mers-based GWAS. The pipeline also contains post-GWAS analysis, such as finding and aligning the source reads for k-mers, aligning kmers and source reads to a reference genome.

## Usage
