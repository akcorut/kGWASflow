# =======================================================================================================
#    Creating a symbolic link for the reference genome fasta
# =======================================================================================================

rule genome_symlink:
    input:
        reference = config["ref"]["fasta"]
    output:
        "resources/ref/genome/genome.fasta",
    message:
        "Creating symbolic link for reference genome fasta..."
    threads: 1
    shell:
        """
        ln -s {input.reference} {output}
        """

# =======================================================================================================
#    Creating bowtie2 index of the reference genome fasta
# =======================================================================================================

rule bowtie2_build:
    input:
        reference = rules.genome_symlink.output,
    output:
        multiext(
            "resources/ref/genome/genome",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 
        config["params"]["align_reads"]["threads"]
    message:
        "Creating bowtie2 index of the reference genome..."
    wrapper:
        "0.80.0/bio/bowtie2/build"

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