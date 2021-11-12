# =======================================================================================================
#     Align reads with significant k-mers to the reference genome
# =======================================================================================================

rule align_reads:
    input:
        reads = rules.fetch_source_reads.output,
        done = "results/fetch_reads_with_kmers/fetch_source_reads.done",
    output:
        dir = directory("results/align_reads_with_kmers/{pheno_filt}"),
        done = "results/align_reads_with_kmers/{pheno_filt}/align_reads.done"
    params:
        in_prefix = lambda w, input: os.path.dirname(input.done),
        out_prefix = lambda w, output: os.path.dirname(output.dir),
        ref_gen = config["ref_gen"]["prefix"],
        pheno = "{pheno_filt}"
    conda:
        "../envs/align_reads.yaml"
    threads: 
        config["params"]["align_reads"]["threads"]
    message:
        ""
    shell:
        """
        python -u scripts/align_reads_with_kmers.py \
        -i {params.in_prefix} \
        -p {params.pheno} \
        -x {params.ref_gen} \
        -o {params.out_prefix} \
        -t {threads}

        touch {output.done}
        """