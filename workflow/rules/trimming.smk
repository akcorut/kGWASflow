# =================================================================================================
#     Trimming for Paired End Reads
# =================================================================================================

rule cutadapt_pe:
    input:
        unpack(get_reads),
    output:
        fastq1="results/trimmed/{sample}/{library}_1.fastq.gz",
        fastq2="results/trimmed/{sample}/{library}_2.fastq.gz",
        qc="results/trimmed/{sample}/{library}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}/{library}.log",
    params:
        others=config["params"]["cutadapt"]["pe"]["extra"],
        adapters=config["params"]["cutadapt"]["pe"]["adapters"],
    threads: 
        config["params"]["cutadapt"]["threads"]
    message: 
        "Performing trimming using cutadapt on on the following files: {input}"
    wrapper:
        "0.79.0/bio/cutadapt/pe"

# =================================================================================================
#     Trimming for Single End Reads
# =================================================================================================

rule cutadapt_se:
    input:
        unpack(get_reads),
    output:
        fastq="results/trimmed/{sample}/{library}.single.fastq.gz",
        qc="results/trimmed/{sample}/{library}.single.qc.txt",
    log:
        "logs/cutadapt/{sample}/{library}.log",
    params:
        others=config["params"]["cutadapt"]["se"]["extra"],
        adapters=config["params"]["cutadapt"]["se"]["adapters"],
    threads: 
        config["params"]["cutadapt"]["threads"]
    message: 
        "Performing trimming using cutadapt on on the following files: {input}"
    wrapper:
        "0.79.0/bio/cutadapt/se"