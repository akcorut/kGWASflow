# =========================================================================================================
#     Filter k-mers (Get the presence/absence patterns for significant k-mers (5per))
# =========================================================================================================

rule filter_kmers:
    input:
        lists = "results/fetch_kmers/{phenos_filt}_kmers_list.txt",
        kmers_table = "results/kmers_table/kmers_table.table",
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        "results/filter_kmers/{phenos_filt}_kmers_table.txt",
    params:
        prefix = lambda w, input: os.path.splitext(input.kmers_table)[0]
    conda:
        "../envs/kmers_gwas.yaml"
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
        
        {input.kmersGWAS_bin}/filter_kmers -t {params.prefix} -k {input.lists} -o {output}
        """
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.fetch_kmers_from_res_table.get(**wildcards).output[0]
    return expand("results/filter_kmers/{phenos_filt}_kmers_table.txt",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

rule check_filter_kmers:
    input:        
        aggregate_input
    output:
        "results/filter_kmers/filter_kmers.done"
    shell:
        """
        touch {output}
        """

