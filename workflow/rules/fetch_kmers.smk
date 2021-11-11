# =========================================================================================================
#     Retrieve the significant k-mers (kmers_pass_threshold_5per) from the results table for each phenotype
# =========================================================================================================

# checkpoint make_dirs_for_fetch_kmers:
#     input:
#         phenos_list= rules.generate_results_table.output.phenos_list_5per
#     output:
#         directory("results/fetch_kmers")
#     params:
#         prefix = "results/fetch_kmers"
#     conda:
#         "../envs/kmers_stats.yaml"
#     shell:
#         """
#         python scripts/make_dirs_for_fetch_kmers.py -p {input.phenos_list} -o {params.prefix}
#         """

checkpoint fetch_kmers_from_res_table:
    input:
        res_tab = rules.generate_results_table.output.res_tab_5per,
        phenos_list= rules.generate_results_table.output.phenos_list_5per,
        # in_dir = rules.make_dirs_for_fetch_kmers.output
    output:
        directory("results/fetch_kmers")
    params:
        # out_prefix= lambda w, output: os.path.dirname(output[0]),
        threshold = config["params"]["fetch_kmers"]["threshold"]
    conda:
        "../envs/kmers_stats.yaml"
    shell:
        """
        python scripts/retrieve_kmers_from_results_table.py -i {input.res_tab} \
        -t {params.threshold} -p {input.phenos_list} -o {output}
        """

# def aggregate_input_fetch_kmers(wildcards):
#     checkpoint_output = checkpoints.make_dirs_for_fetch_kmers.get(**wildcards).output[0]
#     return expand("results/fetch_kmers/{phenos_filt}/{phenos_filt}_kmers_list.txt",
#            phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}")).phenos_filt)


# =========================================================================================================

