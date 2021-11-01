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
#     Min Version and Report
# =================================================================================================

# Make Sure Minimun Snakemake version
min_version("6.10")
basedir = workflow.basedir

# Description of the workflow can be found in the final report
report: "../report/workflow.rst"

# =================================================================================================
#     Laod Config File
# =================================================================================================

configfile: "../config/config.yaml"

# validate(config, schema="../schemas/config.schema.yaml")

# =================================================================================================
#     Sample Sheets and Wildcard Constraints
# =================================================================================================

# samples = (
#     pd.read_csv(
#         config["samples"],
#         sep="\t",
#         dtype={"sample_name": str},
#         comment="#",
#     )
#     .set_index("sample_name", drop=False)
#     .sort_index()
# )

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

# Make sure indes always str type
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
#     Common Helper Functions
# =================================================================================================

def get_fastq(wildcards):
    """Get fastq files using samples sheet."""
    fastqs = samples.loc[(wildcards.sample, wildcards.library), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}

def get_input_path_for_input_list(wildcards):
    if config["settings"]["trimming"]["activate"]:
        return "results/trimmed"
    else:
        return "results/reads"

def get_output_path_for_input_list(wildcards):
    if config["settings"]["trimming"]["activate"]:
        return "results/trimmed/{sample}/input_files.txt"
    else:
        return "results/reads/{sample}/input_files.txt"

