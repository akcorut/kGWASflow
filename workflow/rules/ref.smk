# =======================================================================================================
#    Creating a symbolic link for the reference genome fasta
# =======================================================================================================

rule genome_symlink:
    input:
        reference = config["ref"]["fasta"]
    output:
        "resources/ref/genome/genome.fasta",
    log:
        "logs/ref/genome_symlink/genome_symlink.log"
    message:
        "Creating symbolic link for reference genome fasta..."
    threads: 1
    shell:
        """
        ln -rs {input.reference} {output} 2> {log}
        """

# =======================================================================================================
#    Creating bowtie2 index of the reference genome fasta
# =======================================================================================================

rule bowtie2_build:
    input:
        reference = "resources/ref/genome/genome.fasta",
    output:
        multiext(
            "resources/ref/genome/bowtie2_index/genome",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/ref/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 
        config["params"]["bowtie2"]["threads"]
    message:
        "Creating bowtie2 index of the reference genome..."
    wrapper:
        "0.80.0/bio/bowtie2/build"

# =======================================================================================================
#    Creating bowtie index of the reference genome fasta
# =======================================================================================================

rule bowtie_build:
    input:
        reference = "resources/ref/genome/genome.fasta",
    output:
        multiext(
            "resources/ref/genome/bowtie_index/genome",
            ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt",
        ),
    log:
        "logs/bowtie_build/build.log"
    params:
        prefix= lambda w, output: os.path.commonprefix(output).rstrip("."),
        extra=""  # optional parameters
    conda:
        "../envs/align_kmers.yaml"
    threads: 
        config["params"]["bowtie"]["threads"]
    log:
        "logs/ref/bowtie_build/build.log"
    message:
        "Creating bowtie index of the reference genome..."
    shell:
        "bowtie-build --threads {threads} {input.reference} {params.prefix} 2> {log}"

# =======================================================================================================
#    SRA-download
# =======================================================================================================

rule sra_get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "resources/ref/sra-pe-reads/{accession}_1.fastq",
        "resources/ref/sra-pe-reads/{accession}_2.fastq"
    params:
        extra="--skip-technical --progress --temp resources/ref/sra-pe-reads"
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "v1.12.2/bio/sra-tools/fasterq-dump"

rule sra_get_fastq_se:
    output:
        "resources/ref/sra-se-reads/{accession}.fastq"
    params:
        extra="--skip-technical --progress --temp resources/ref/sra-pe-reads"
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "v1.12.2/bio/sra-tools/fasterq-dump"