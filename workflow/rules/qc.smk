# =================================================================================================
#     Generate QC Stats Using FastQC on Raw Reads
# =================================================================================================

rule fastqc:
    input:
        unpack(get_reads)
    output:
        html="results/qc/fastqc/{sample}/{library}_fastqc.html",
        zip="results/qc/fastqc/{sample}/{library}_fastqc.zip",
    log:
        "logs/fastqc/{sample}/{library}.log",
    threads:
        config["params"]["fastqc"]["threads"]
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
            "results/qc/fastqc/{u.sample_name}/{u.library_name}_fastqc.zip",
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
    params:
        "-v -d"
    message: 
        "Performing MultiQC on the FastQC results..."
    wrapper:
        "0.79.0/bio/multiqc"

# =================================================================================================