rule bowtie2_build:
    input:
        reference= config["ref"]["fasta"]
    output:
        multiext(
            lambda w, input: os.path.splitext(input.reference)[0],
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "0.80.0/bio/bowtie2/build"