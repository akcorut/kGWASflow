# =======================================================================================================
#    Creating a symbolic link for the reference genome fasta
# =======================================================================================================

rule genome_symlink:
    input:
        reference = config["ref"]["fasta"]
    output:
        "resources/genome.fasta",
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
            "resources/genome",
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