# =========================================================================================================
#     Retrieve the significant k-mers (kmers_pass_threshold_5per) from the results table for each phenotype
# =========================================================================================================

checkpoint fetch_kmers_from_res_table:
    input:
        res_tab = rules.generate_results_table.output.res_tab_5per,
        phenos_list= rules.generate_results_table.output.phenos_list_5per,
    output:
        directory("results/fetch_kmers")
    params:
        threshold = config["params"]["fetch_kmers"]["threshold"]
    conda:
        "../envs/kmers_stats.yaml"
    message:
        "Fetching significant k-mers for each phenotype from kmersGWAS results table ({input.res_tab})..."
    shell:
        """
        python scripts/fetch_kmers_from_results_table.py -i {input.res_tab} \
        -t {params.threshold} -p {input.phenos_list} -o {output}
        """

# =========================================================================================================

