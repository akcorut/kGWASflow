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
    threads: 8
    wrapper:
        "0.80.0/bio/bowtie2/build"