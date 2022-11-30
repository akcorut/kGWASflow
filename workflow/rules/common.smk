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
from datetime import datetime

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
sample_names=samples.index.get_level_values("sample_name")
library_names=samples.index.get_level_values("library_name")

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

kgwasflow_version = "v1.0.0"
kgwasflow_author = "Adnan Kivanc Corut"
date_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
snake_version = snakemake.__version__
py_version= sys.version.split(' ')[0]
snakefile_path= workflow.snakefile
base_dir= workflow.basedir
work_dir= os.getcwd()
config_file= ", ".join(workflow.configfiles)

# Helpful messages
logger.info("# ================================================================================== #")
logger.info("      _     _______          __      _____  __ _                 ") 
logger.info("     | |   / ____\ \        / /\    / ____|/ _| |                ")
logger.info("     | | _| |  __ \ \  /\  / /  \  | (___ | |_| | _____      __  ")
logger.info("     | |/ / | |_ | \ \/  \/ / /\ \  \___ \|  _| |/ _ \ \ /\ / /  ")
logger.info("     |   <| |__| |  \  /\  / ____ \ ____) | | | | (_) \ V  V /   ")
logger.info("     |_|\_\\_____|   \/  \/_/    \_\_____/|_| |_|\___/ \_/\_/    ")                                      
logger.info("")
logger.info("     kGWASflow: A Snakemake Workflow for k-mers Based GWAS                            ")
logger.info("                                                                                      ")
logger.info(f"     Date:            {date_time}")
logger.info(f"     Author:          {kgwasflow_author}")
logger.info(f"     kGWASflow version:          {kgwasflow_version}")
logger.info(f"     Snakemake version:          {snake_version}")
logger.info(f"     Python version:             {py_version}")
logger.info("")
logger.info(f"     Snakefile:          {snakefile_path}")
logger.info(f"     Base directory:     {base_dir}")
logger.info(f"     Working directory:  {work_dir}")
logger.info(f"     Config files:       {config_file}")
logger.info("                                                                                      ")
logger.info("# ================================================================================== #")
logger.info("")
logger.info("     If you use kGWASflow, please cite:                                            ")
logger.info("")
logger.info("     * Voichek, Y., Weigel, D. Identifying genetic variants underlying                  ")
logger.info("     phenotypic variation in plants without complete genomes.                         ")
logger.info("     Nat Genet 52, 534–540 (2020). https://doi.org/10.1038/s41588-020-0612-7          ")
logger.info("")
logger.info("     * Kivanc Corut. akcorut/kGWASflow: v1.0.0. (2022). ")
logger.info("     https://doi.org/10.5281/zenodo.7290926                                        ")
logger.info("")
logger.info("# ================================================================================== #")
logger.info("")

# =================================================================================================
#     Set kmersGWAS Version and Path
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

def sra_only(sample, library):
    return pd.isnull(samples.loc[(sample, library), "fq1"]) and pd.isnull(samples.loc[(sample, library), "fq2"]) \
           and not pd.isnull(samples.loc[(sample, library), "sra"])

def is_sra_se(sample, library):
    return sra_only(sample, library) and config["settings"]["single_end"]

def is_sra_pe(sample, library):
    return sra_only(sample, library) and not config["settings"]["single_end"]

def get_individual_fastq(wildcards):
# Adapted from: https://github.com/snakemake-workflows/chipseq
    """Get individual raw FASTQ files from samples sheet."""
    fastqs = samples.loc[(wildcards.sample, wildcards.library), ["fq1", "fq2"]].dropna()
    if ( len(fastqs) == 0 or len(fastqs) == 1 ):
        if is_sra_se(wildcards.sample, wildcards.library):
            return expand("resources/ref/sra-se-reads/{accession}.fastq",
                              accession=samples.loc[ (wildcards.sample, wildcards.library), "sra" ])
        elif is_sra_pe(wildcards.sample, wildcards.library):
            return expand("resources/ref/sra-pe-reads/{accession}_1.fastq",
                              accession=samples.loc[ (wildcards.sample, wildcards.library), "sra" ])
        else:
            return samples.loc[ (wildcards.sample, wildcards.library), "fq1" ]
    elif len(fastqs) == 2:
        if is_sra_pe(wildcards.sample, wildcards.library):
            return expand("resources/ref/sra-pe-reads/{accession}_2.fastq",
                          accession=samples.loc[ (wildcards.sample, wildcards.library), "sra" ])
        else:
            return samples.loc[ (wildcards.sample, wildcards.library), "fq2" ]

def get_fastqs(wildcards):
# Adapted from: https://github.com/snakemake-workflows/chipseq
    """Get raw FASTQ files from samples sheet."""
    fastqs = samples.loc[(wildcards.sample, wildcards.library), ["fq1", "fq2"]].dropna()
    if is_sra_se(wildcards.sample, wildcards.library):
        return expand("resources/ref/sra-se-reads/{accession}.fastq",
                          accession=samples.loc[ (wildcards.sample, wildcards.library), "sra" ])
    elif is_sra_pe(wildcards.sample, wildcards.library):
        return expand(["resources/ref/sra-pe-reads/{accession}_1.fastq", "resources/ref/sra-pe-reads/{accession}_2.fastq"],
                          accession=samples.loc[ (wildcards.sample, wildcards.library), "sra" ])
    elif len(fastqs) == 1:
        return samples.loc[ (wildcards.sample, wildcards.library), "fq1" ]
    else:
        u = samples.loc[ (wildcards.sample, wildcards.library), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def ends_with_gz(samp_tab):
    fqs = samp_tab.loc[:, ['fq1']]
    if not pd.isnull(fqs.fq1).all():
        if samp_tab["fq1"].str.endswith('gz').all():
            return True
        else:
            return False
    else:
        return False

def get_multiqc_input(wildcards):
    """Get multiqc input."""
    multiqc_input = []
    multiqc_input.extend(
        expand(
            [
                "results/qc/fastqc/{u.sample_name}/{u.library_name}_fastqc.zip",
            ],
            u=samples.itertuples()
        )
    )
    if config["settings"]["trimming"]["activate"]:
        multiqc_input.extend(
            expand(
                [
                    "results/trimmed/{u.sample_name}/{u.library_name}.paired.qc.txt",
                ],
                u=samples.itertuples()
            )
        )
    return multiqc_input

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

def get_plink_prefix():
    """Get prefix fopr the SNPs plink file."""
    plink_path = config["settings"]["kmers_gwas"]["use_snps_kinship"]["snps_plink"]
    plink_prefix = os.path.splitext(plink_path)[0]
    return plink_prefix

# ================= Checkpoint Fucntions ================= #

def aggregate_input_filter_kmers(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/filter_kmers/{phenos_filt}_kmers_table.txt",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_fetch_reads(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/fetch_reads_with_kmers/{phenos_filt}",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_contigs(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam.bai",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_kmers(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    if not config["settings"]["align_kmers"]["plot_manhattan"]:
        return expand("results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam.bai",
               phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)
    else:
        return expand(
                   [
                        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam.bai",
                        "results/plots/manhattan/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.manhattan_plot.pdf",
                   ], phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt
                    )

def aggregate_input_align_reads(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam.bai",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_blast_contigs(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/blast_contigs/{phenos_filt}/{phenos_filt}_contigs.blast.txt",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

# =================================================================================================
#     Target Ouput Function
# =================================================================================================

def get_target_output(wildcards):
    """
    Get all requested inputs (target outputs) for rule all.
    """

    align_kmers = config["settings"]["align_kmers"]["activate"]
    align_reads_with_kmers = config["settings"]["align_reads_with_kmers"]["activate"]
    assemble_reads_with_kmers = config["settings"]["assemble_reads"]["activate"]
    blast_assembled_reads = config["settings"]["blast_contigs"]["activate"]

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

    if align_kmers:
        target_output.extend(
            expand(
                "results/align_kmers/align_kmers.done"
            )
        ),
    
    if align_reads_with_kmers:
        target_output.extend(
            expand(
                "results/align_reads_with_kmers/align_reads_with_kmers.done"
            )
        ),

    if assemble_reads_with_kmers:
        target_output.extend(
            expand(
                "results/align_contigs/align_contigs.done"
            )
        ),

    if blast_assembled_reads:
        target_output.extend(
            expand(
                "results/blast_contigs/blast_contigs.done"
            )
        ),

    return target_output

# =================================================================================================