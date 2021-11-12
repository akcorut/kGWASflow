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

def aggregate_input_align_reads(wildcards):
    checkpoint_output = checkpoints.fetch_kmers_from_res_table.get(**wildcards).output[0]
    return expand("results/align_reads_with_kmers/{pheno_filt}",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

rule check_fetch_reads:
    input:
        aggregate_input_align_reads
    output:
        "results/align_reads_with_kmers/align_reads_with_kmers.done"
    message:
        "Checking if aligning reads with significant k-mers is done..."
    shell:
        """
        touch {output}
        """