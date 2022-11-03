# =========================================================================================================
#     Filter k-mers (Get the presence/absence patterns for significant k-mers (5per))
#     # Source code of 'filter_kmers': https://github.com/voichek/kmersGWAS
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
    message:
        "Filtering k-mers table using {input.lists} and output the results in textual format..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib

        {input.kmersGWAS_bin}/filter_kmers -t {params.prefix} -k {input.lists} -o {output}
        """

# =========================================================================================================
#     Aggregate filter_kmers outputs 
# =========================================================================================================

rule aggregate_filter_kmers:
    input:        
        aggregate_input_filter_kmers
    output:
        "results/filter_kmers/filter_kmers.done"
    message:
        "Checking if filtering k-mers steps is done..."
    shell:
        """
        touch {output}
        """

# =========================================================================================================
