# =======================================================================================================
#     Download and compile fetch_reads_with_kmers 
#     # Source: https://github.com/voichek/fetch_reads_with_kmers
# =======================================================================================================


rule download_fetch_reads_with_kmers:
    output:
        cpp= "scripts/external/fetch_reads_with_kmers/fetch_reads.cpp",
    shell:
        """
        git clone --recurse-submodules https://github.com/voichek/fetch_reads_with_kmers.git
        mv fetch_reads_with_kmers scripts/external
        """

rule make_fetch_reads_with_kmers:
    input:
        rules.download_fetch_reads_with_kmers.output.cpp
    output:
        "scripts/external/fetch_reads_with_kmers/fetch_reads"
    params:
        prefix= lambda w, output: os.path.dirname(output[0])
    shell:
        """
        cd {params.prefix}
        make
        """

# =======================================================================================================
#     Fetch source reads of significant k-mers
# =======================================================================================================

rule fetch_source_reads:
    input:
        kmers_tab = "results/filter_kmers/{phenos_filt}_kmers_table.txt",
        kmers_list = "results/fetch_kmers/{phenos_filt}_kmers_list.fa",
        fetch_reads = "scripts/external/fetch_reads_with_kmers/fetch_reads",
        samp_tab = config["samples"],
        filter_kmers = "results/filter_kmers/filter_kmers.done"
    output:
        directory("results/fetch_reads_with_kmers/{phenos_filt}"),
    params:
        kmers_tab_prefix = lambda w, input: os.path.dirname(input.kmers_tab),
        kmers_list_prefix = lambda w, input: os.path.dirname(input.kmers_list),
        pheno = "{phenos_filt}",
        out_prefix = "results/fetch_reads_with_kmers",
        fetch_reads_prefix = lambda w, input: os.path.dirname(input.fetch_reads)
    conda:
        "../envs/kmers_stats.yaml"
    log:
        "logs/fetch_reads_with_kmers/{phenos_filt}/fetch_source_reads_of_kmers.log"
    message:
        "Fetching reads that contain significant k-mers find in {input.kmers_list}..."    
    script:
        "../scripts/fetch_source_reads.py"

# =========================================================================================================
#     Check fetch_source_reads 
# =========================================================================================================

def aggregate_input_fetch_reads(wildcards):
    checkpoint_output = checkpoints.fetch_kmers_from_res_table.get(**wildcards).output[0]
    return expand("results/fetch_reads_with_kmers/{phenos_filt}",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

rule check_fetch_reads:
    input:
        aggregate_input_fetch_reads
    output:
        "results/fetch_reads_with_kmers/fetch_source_reads.done"
    message:
        "Checking if fetching source reads of k-mers is done..."
    shell:
        """
        touch {output}
        """

# =========================================================================================================