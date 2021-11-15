# =======================================================================================================
#     Align reads with significant k-mers to the reference genome
# =======================================================================================================

rule align_reads:
    input:
        index = rules.bowtie2_build.output,
        reads = rules.fetch_source_reads.output,
        done = "results/fetch_reads_with_kmers/fetch_source_reads.done",
    output:
        dir = directory("results/align_reads_with_kmers/{phenos_filt}"),
        done = "results/align_reads_with_kmers/{phenos_filt}/align_reads.done"
    params:
        in_prefix = lambda w, input: os.path.dirname(input.done),
        out_prefix = lambda w, output: os.path.dirname(output.dir),
        pheno = "{phenos_filt}",
        index = "resources/genome",
    conda:
        "../envs/align_reads.yaml"
    threads: 
        config["params"]["align_reads"]["threads"]
    message:
        "Aligning reads with k-mers to the reference genome..."
    shell:
        """
        python -u scripts/align_reads_with_kmers.py \
        -i {params.in_prefix} \
        -p {params.pheno} \
        -r {params.index} \
        -o {params.out_prefix} \
        -t {threads}

        touch {output.done}
        """

# =========================================================================================================
#     Check align_reads 
# =========================================================================================================

def aggregate_input_align_reads(wildcards):
    checkpoint_output = checkpoints.fetch_kmers_from_res_table.get(**wildcards).output[0]
    return expand("results/align_reads_with_kmers/{phenos_filt}",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

rule check_align_reads:
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

# =========================================================================================================