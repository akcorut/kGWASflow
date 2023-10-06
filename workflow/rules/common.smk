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

kgwasflow_version = "v1.3"
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
logger.info("     * Corut, A. K. & Wallace, J. G. kGWASflow: a modular, flexible,                   ")
logger.info("     and reproducible Snakemake workflow for k-mers-based GWAS.                        ")
logger.info("     bioRxiv (2023). https://doi.org/10.1101/2023.07.10.548365                         ")
logger.info("")
logger.info("     * Voichek, Y., Weigel, D. Identifying genetic variants underlying                  ")
logger.info("     phenotypic variation in plants without complete genomes.                         ")
logger.info("     Nat Genet 52, 534–540 (2020). https://doi.org/10.1038/s41588-020-0612-7          ")
logger.info("")
logger.info("# ================================================================================== #")
logger.info("")

# =================================================================================================
#     Set kmersGWAS Version and Path
#     # Source: https://github.com/voichek/kmersGWAS/
#     # Reference: https://www.nature.com/articles/s41588-020-0612-7
# =================================================================================================

KMERSGWAS_DIR = "scripts/external/kmers_gwas"
KMERSGWAS_VERSION = "v0.3-beta"
KMERSGWAS_ZIP_PREFIX = "v0_3_beta"
KMERSGWAS_ZIP_NAME = f"{KMERSGWAS_VERSION}.zip"
KMERSGWAS_ZIP_PATH = os.path.join(KMERSGWAS_DIR, KMERSGWAS_ZIP_NAME)
KMERSGWAS_PY_PATH = os.path.join(KMERSGWAS_DIR, "kmers_gwas.py")
KMERSGWAS_BIN_PATH = os.path.join(KMERSGWAS_DIR, "bin/")

# =================================================================================================
#     Common Helper Functions
# =================================================================================================

# Check if the sample is SRA only
def sra_only(sample, library): 
    return pd.isnull(samples.loc[(sample, library), "fq1"]) and pd.isnull(samples.loc[(sample, library), "fq2"]) \
           and not pd.isnull(samples.loc[(sample, library), "sra"])

# Check if the sample is SRA only and single-end
def is_sra_se(sample, library): 
    return sra_only(sample, library) and config["settings"]["single_end"]

# Check if the sample is SRA only and paired-end
def is_sra_pe(sample, library): 
    return sra_only(sample, library) and not config["settings"]["single_end"]

# Get individual FASTQ files
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

# Get FASTQ files
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

# Check if FASTQ files are gzipped
def ends_with_gz(samp_tab):
    fqs = samp_tab.loc[:, ['fq1']]
    if not pd.isnull(fqs.fq1).all():
        if samp_tab["fq1"].str.endswith('gz').all():
            return True
        else:
            return False
    else:
        return False

# Get MultiQC input
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

# Get phenotypes
def get_phenos(wildcards): 
    """Get fastq files using samples sheet."""
    pheno_names = phenos.loc[(wildcards.pheno), ["pheno_path"]].dropna()
    return {"pheno_path": pheno_names.pheno_path}

def get_plink_prefix():
    """Get prefix fopr the SNPs plink file."""
    # if plink file is not provided exit and print error message
    if config["settings"]["kmers_gwas"]["use_snps_kinship"]["snps_plink"] == "":
        print("Please provide a PLINK file in the config file.")
        sys.exit(1)
    else:
        plink_path = config["settings"]["kmers_gwas"]["use_snps_kinship"]["snps_plink"]
        plink_prefix = os.path.splitext(plink_path)[0]
        return plink_prefix

def get_referece_genome():
    """Get reference genome."""
    # if reference genome is not provided exit and print error message
    if config["ref"]["fasta"] == "":
        print("Please provide a reference genome FASTA file in the config file.")
        sys.exit(1)
    else:
        return config["ref"]["fasta"]

def get_genome_annotation():
    """Get genome annotation."""
    # if genome annotation is not provided exit and print error message
    if config["ref"]["annotation"] == "":
        print("Please provide a genome annotation GTF/GFF file in the config file.")
        sys.exit(1)
    else:
        return config["ref"]["annotation"]

# ================= Checkpoint Fucntions ================= #

def aggregate_input_filter_kmers(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/filter_kmers/{phenos_filt}_kmers_table.txt",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_fetch_reads(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/fetch_reads_with_kmers/{phenos_filt}/individual_reads",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_kmers(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand([
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam.bai",
        "results/plots/manhattan/align_kmers/{phenos_filt}/{phenos_filt}.kmers_alignment.manhattan_plot.pdf",
    ], phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_reads(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.merged.bed",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_contigs(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand([
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bed",
    ], phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_blast_contigs(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/blast_contigs/{phenos_filt}/{phenos_filt}_contigs.blast.txt",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_contigs_igv_report(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/igv_reports/align_contigs/{phenos_filt}/{phenos_filt}_contigs_aligned.igv-report.html",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

def aggregate_input_align_reads_igv_report(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/igv_reports/align_reads/{phenos_filt}/{phenos_filt}_reads_aligned.igv-report.html",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

# =================================================================================================
#     Target Ouput Function
# =================================================================================================

def get_target_output(wildcards):
    """
    Get all requested inputs (target outputs) for rule all.
    """

    # Activate/Deactivate rules
    align_kmers = config["settings"]["align_kmers"]["activate"] # Activate align kmers
    align_reads_with_kmers = config["settings"]["align_reads_with_kmers"]["activate"] # Activate align reads with kmers
    assemble_reads_with_kmers = config["settings"]["assemble_reads"]["activate"] # Activate assemble reads with kmers
    blast_assembled_reads = config["settings"]["blast_contigs"]["activate"] # Activate BLAST assembled reads
    convert_kmers_table_to_plink = config["settings"]["kmers_gwas"]["kmers_table_to_bed"]["activate"] # Activate convert kmers table to plink
    align_contigs_igv_report = config["settings"]["assemble_reads"]["igv_report"] # Activate align contigs IGV report
    align_reads_igv_report = config["settings"]["align_reads_with_kmers"]["igv_report"] # Activate align reads IGV report

    target_output = [] # List of target outputs

    ### Add target outputs for rule all ###
    
    # MultiQC
    target_output.extend(
        expand(
            "results/qc/multiqc.html"
        )
    ),

    # KMC, k-mers stats
    target_output.extend(
        expand(
            [
                "results/plots/kmer_stats/kmc_canon.total_reads_vs_unique_kmers.joint_plot.pdf",
                "results/plots/kmer_stats/kmc_all.total_reads_vs_unique_kmers.joint_plot.pdf",
                "results/plots/kmer_stats/kmc_canon.total_reads_vs_unique_kmers.scatter_plot.pdf",
                "results/plots/kmer_stats/kmc_all.total_reads_vs_unique_kmers.scatter_plot.pdf",
                "results/plots/kmer_dist/kmc_canon.kmer_dist_hist.pdf",
                "results/plots/kmer_dist/kmc_all.kmer_dist_hist.pdf",
                "results/tables/kmer_stats/kmc_canon.kmer_stats.tsv",
                "results/tables/kmer_stats/kmc_all.kmer_stats.tsv"
            ]
        )
    ),

    # Combine and filter k-mers
    target_output.extend(
        expand(
            "results/plots/kmer_allele_counts/kmer_allele_counts.barplot.pdf",
        )
    ),

    #  Convert the k-mers table to PLINK
    if convert_kmers_table_to_plink: 
        # If convert_kmers_table_to_plink is activated
        target_output.extend(
            expand(
               "results/kmers_table/plink/pheno_{pheno}/convert_kmers_table_to_plink.done",
                pheno= config["params"]["kmers_table_to_bed"]["phenos"]
            )
        ),

    # kmersGWAS
    target_output.extend(
        expand(
            "results/kmers_gwas/{pheno}/kmers/pass_threshold_5per",
            pheno= pheno_names
        )
    ),

    # kmersGWAS results summary
    target_output.extend(
        expand(
            "results/tables/kmers_gwas_summary/{pheno}/generate_results_summary.done",
            pheno= pheno_names
        )
    ),

    # Align k-mers
    if align_kmers:
        target_output.extend(
            expand(
                "results/align_kmers/align_kmers.done"
            )
        ),
        
    # Align reads with k-mers
    if align_reads_with_kmers:
        target_output.extend(
            expand(
                "results/align_reads_with_kmers/align_reads_with_kmers.done"
            )
        ),

    # Assemble reads with k-mers and align contigs
    if assemble_reads_with_kmers:
        target_output.extend(
            expand(
                "results/align_contigs/align_contigs.done"
            )
        ),
    
    # IGV report for aligned contigs
    if assemble_reads_with_kmers and align_contigs_igv_report:
        target_output.extend(
            expand(
                "results/igv_reports/align_contigs.igv_report.done"
            )
        ),
    
    # IGV report for aligned reads
    if align_reads_with_kmers and align_reads_igv_report:
        target_output.extend(
            expand(
                "results/igv_reports/align_reads.igv_report.done"
            )
        ),

    # BLAST assembled reads
    if blast_assembled_reads:
        target_output.extend(
            expand(
                "results/blast_contigs/blast_contigs.done"
            )
        ),

    # Return target outputs
    return target_output

# =================================================================================================