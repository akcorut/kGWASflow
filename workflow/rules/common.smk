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
#     Laod Config File
# =================================================================================================

# validate(config, schema="../schemas/config.schema.yaml")

# =================================================================================================
#     Sample Sheets and Wildcard Constraints
# =================================================================================================

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

# Make sure indeces always str type
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])

# Generate sample and library names lists
sample_names=list(set(samples.index.get_level_values("sample_name")))
library_names=list(set(samples.index.get_level_values("library_name")))

# Wildcard constraints: Allowing only sample and library names from the samples sheet to be used
wildcard_constraints:
    sample="|".join(sample_names),
    library="|".join(library_names)

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

## Source: https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/40e9aed1e2534ae6d262b6f7a05e8b625f127682/workflow/rules/common.smk#L59-L62
def is_single_end(sample, library):
    no_fq2 = pd.isnull(samples.loc[(sample, library), "fq2"])
    return no_fq2

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
            "results/kmers_count/{sample}/kmers_with_strand", sample=sample_names
        )
    ),
    target_output.extend(
        expand(
            [
                "results/plots/kmers_count/kmc_canon_total_reads_vs_unique_kmers.joint_plot.png",
                "results/plots/kmers_count/kmc_all_total_reads_vs_unique_kmers.joint_plot.png"
            ]
        )
    )
    return target_output