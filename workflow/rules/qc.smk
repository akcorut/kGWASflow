# =================================================================================================
#     Generate QC Stats Using FastQC on Raw Reads
# =================================================================================================

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="results/qc/fastqc/{sample}/{library}.html",
        zip="results/qc/fastqc/{sample}/{library}.zip",
    log:
        "logs/fastqc/{sample}/{library}.log",
    message: 
        "Performing quality control analysis using FastQC on the following files: {input}"
    wrapper:
        "0.79.0/bio/fastqc"

# =================================================================================================
#     MultiQC
# =================================================================================================

rule multiqc:
    input:
        expand(
            "results/qc/fastqc/{u.sample_name}/{u.library_name}.zip",
            u=samples.itertuples(),
        )
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    message: 
        "Performing MultiQC on the FastQC results..."
    wrapper:
        "0.79.0/bio/multiqc"

# =================================================================================================
#     FastQC After Trimming
# =================================================================================================



# =================================================================================================
#     MultiQC After Trimming
# =================================================================================================