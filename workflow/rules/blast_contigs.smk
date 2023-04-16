# =======================================================================================================
#     Create a BLAST database
# =======================================================================================================

rule blast_makedb_nuc:
    input:
        fasta= "resources/ref/genome/genome.fasta"
    output:
        multiext("resources/ref/genome/genome.fasta",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    log:
        "logs/blast_contigs/makeblastdb.log"
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids"
    wrapper:
        "v1.25.0/bio/blast/makeblastdb"

# =======================================================================================================
#     Filter out short contigs
# =======================================================================================================

rule filter_contigs:
    input:
        "results/assemble_reads_with_kmers/{phenos_filt}/assembly/contigs.fasta"
    output:
        "results/assemble_reads_with_kmers/{phenos_filt}/filtered_contigs/{phenos_filt}_contigs.filtered.fasta"
    conda:
        "../envs/align_reads.yaml"
    log:
        "logs/blast_contigs/filter_contigs/filter.{phenos_filt}_contigs.log"
    shell:
        """
        seqkit head -n 1 {input} > {output} 2> {log}
        """

# =======================================================================================================
#     BLAST longest contigs against reference db
# =======================================================================================================

rule blast_contigs:
    input:
        query = "results/assemble_reads_with_kmers/{phenos_filt}/filtered_contigs/{phenos_filt}_contigs.filtered.fasta",
        blastdb=multiext("resources/ref/genome/genome.fasta",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    output:
        "results/blast_contigs/{phenos_filt}/{phenos_filt}_contigs.blast.txt"
    log:
        "logs/blast_contigs/{phenos_filt}_contigs.blast.log"
    threads: 
        config["params"]["blastn"]["threads"]
    params:
        format= config["params"]["blastn"]["format"],
        extra= config["params"]["blastn"]["extra"]
    wrapper:
        "v1.25.0/bio/blast/blastn"

# =========================================================================================================
#     Aggregate blast_contigs outputs
# =========================================================================================================

rule aggregate_blast_contigs:
    input:
        aggregate_input_blast_contigs
    output:
        "results/blast_contigs/blast_contigs.done"
    log:
        "logs/blast_contigs/aggregate_blast_contigs.log"
    shell:
        """
        touch {output}
        """

# =========================================================================================================