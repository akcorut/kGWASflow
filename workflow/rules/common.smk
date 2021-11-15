# =================================================================================================
#     Import Libraries
# =================================================================================================

import glob
import os
import json
import pandas as pd
import socket, platform
from snakemake.io import expand
from snakemake.utils import R
from snakemake.utils import min_version
from snakemake.utils import validate
from snakemake.io import glob_wildcards
import re
from os.path import join, basename, dirname
import pathlib
from os import path

# =================================================================================================
#     Sample and Phenotype Sheets + Wildcard Constraints
# =================================================================================================

# ===================== Samples ===================== #

# Read samples sheet (samples.tsv)
samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "library_name": str},
        comment="#",
    )
    .set_index(["sample_name", "library_name"], drop=False)
    .sort_index()
)
samples.index.names = ["sample_name", "library_name"]

# Validating the samples sheet
validate(samples, schema="../schemas/samples.schema.yaml")

# Make sure indeces always str type
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])

# Generate sample and library names lists
sample_names=list(set(samples.index.get_level_values("sample_name")))
library_names=list(set(samples.index.get_level_values("library_name")))

# ==================== Phenotypes ===================== #

# Read phenotypes sheet (phenos.tsv)
phenos = (
    pd.read_csv(
        config["phenotypes"],
        sep="\t",
        dtype={"pheno_name": str},
        comment="#",
    )
    .set_index(["pheno_name"], drop=False)
    .sort_index()
)

# Validating the phenotpyes sheet
validate(phenos, schema="../schemas/phenos.schema.yaml")

# Generate phenotype names list
pheno_names=list(set(phenos.index.get_level_values("pheno_name")))

# ==================== Wildcards ===================== #

# Wildcard constraints: Allowing only sample and 
# library names from the samples sheet to be used
wildcard_constraints:
    sample="|".join(sample_names),
    library="|".join(library_names),
    pheno="|".join(pheno_names)

# =================================================================================================
#     Pipeline User Output
# =================================================================================================

# Helpful messages
logger.info("===========================================================================")
logger.info("")
logger.info("    kGWASflow: A Snakemake Workflow for k-mers Based GWAS          ")
logger.info("")
logger.info("    Snakefile:          " + (workflow.snakefile))
logger.info("    Base directory:     " + (workflow.basedir))
logger.info("    Config files:       " + (", ".join(workflow.configfiles)))
logger.info("")
logger.info("===========================================================================")
logger.info("")

# =================================================================================================
#     Determine kmersGWAS Version and Path
#     # Source: https://github.com/voichek/kmersGWAS/
#     # Reference: https://www.nature.com/articles/s41588-020-0612-7
# =================================================================================================

KMERSGWAS_DIR = "scripts/external/kmers_gwas"
KMERSGWAS_VERSION = "v0.2-beta"
KMERSGWAS_ZIP_PREFIX = "v0_2_beta"
KMERSGWAS_ZIP_NAME = f"{KMERSGWAS_VERSION}.zip"
KMERSGWAS_ZIP_PATH = os.path.join(KMERSGWAS_DIR, KMERSGWAS_ZIP_NAME)
KMERSGWAS_PY_PATH = os.path.join(KMERSGWAS_DIR, "kmers_gwas.py")
KMERSGWAS_BIN_PATH = os.path.join(KMERSGWAS_DIR, "bin/")

# =================================================================================================
#     Common Helper Functions
# =================================================================================================

def get_reads(wildcards):
    """Get fastq files using samples sheet."""
    fastqs = samples.loc[(wildcards.sample, wildcards.library), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        # paired-end
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        # single end
        return {"r1": fastqs.fq1}

def get_phenos(wildcards):
    """Get fastq files using samples sheet."""
    pheno_names = phenos.loc[(wildcards.pheno), ["pheno_path"]].dropna()
    return {"pheno_path": pheno_names.pheno_path}

def get_input_path_for_generate_input_lists():
    """Get input path of reads to create input lists."""
    if config["settings"]["trimming"]["activate"]:
        return "results/trimmed"
    else:
        return "results/reads"

def get_generate_input_lists_target():
    if config["settings"]["trimming"]["activate"]:
        return expand(
            "results/trimmed/{sample}/input_files.txt", sample=sample_names
        )
    else:
        return expand(
            "results/reads/{sample}/input_files.txt", sample=sample_names
        )

def get_plink_prefix():
    plink_path = config["settings"]["kmers_gwas"]["use_snps_kinship"]["snps_plink_file"]
    plink_prefix = os.path.splitext(plink_path)[0]
    return plink_prefix


# =================================================================================================
#     Target Ouput Function
# =================================================================================================

def get_target_output(wildcards):
    """
    Get all requested inputs (target outputs) for rule all.
    """

    target_output = []

    target_output.extend(
        expand(
            "results/qc/multiqc.html"
        )
    ),
    target_output.extend(
        expand(
            [
                "results/plots/kmers_count/kmc_canon_total_reads_vs_unique_kmers.joint_plot.pdf",
                "results/plots/kmers_count/kmc_all_total_reads_vs_unique_kmers.joint_plot.pdf",
                "results/plots/kmers_count/kmc_canon_total_reads_vs_unique_kmers.scatter_plot.pdf",
                "results/plots/kmers_count/kmc_all_total_reads_vs_unique_kmers.scatter_plot.pdf",
                "results/tables/kmers_count/kmc_canon.stats.tsv",
                "results/tables/kmers_count/kmc_all.stats.tsv"
            ]
        )
    ),
    target_output.extend(
        expand(
            "results/plots/kmers_list/kmer_allele_counts.pdf"
        )
    ),
    target_output.extend(
        expand(
            "results/kmers_gwas/{pheno}/kmers/pass_threshold_5per",
            pheno= pheno_names
        )
    ),
    # target_output.extend(
    #     expand(
    #         "results/fetch_reads_with_kmers/fetch_source_reads.done"
    #     )
    # ),
    target_output.extend(
        expand(
            "results/align_reads_with_kmers/align_reads_with_kmers.done"
        )
    ),
    return target_output

# =================================================================================================