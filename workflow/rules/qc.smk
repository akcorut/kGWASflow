# =================================================================================================
#     Generate QC Stats Using FastQC on Raw Reads
# =================================================================================================

rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}/{library}_fastqc.html",
        zip="results/qc/fastqc/{sample}/{library}_fastqc.zip",
    log:
        "logs/qc/fastqc/{sample}/{library}.log",
    threads:
        config["params"]["fastqc"]["threads"]
    message: 
        "Performing quality control analysis using FastQC on the following files: {input}"
    wrapper:
        "v1.12.2/bio/fastqc"

# =================================================================================================
#     MultiQC
# =================================================================================================

rule multiqc:
    input:
        get_multiqc_input
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/qc/multiqc/multiqc.log",
    params:
        "-v -d --interactive"
    message: 
        "Performing MultiQC on the FastQC results..."
    wrapper:
        "v1.12.2/bio/multiqc"

# =================================================================================================