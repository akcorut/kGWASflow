# =========================================================================================================
#     Retrieve the significant k-mers (kmers_pass_threshold_5per) from the results table for each phenotype
# =========================================================================================================

checkpoint fetch_significant_kmers:
    input:
        res_tab = rules.generate_results_table.output.res_tab_5per,
    output:
        directory("results/fetch_kmers/")
    params:
        in_prefix = "results/kmers_gwas",
        threshold = config["params"]["fetch_kmers"]["threshold"],
        check_file = "results/kmers_gwas_summary/NO_KMERS_PASS_5PER_THRESHOLD_FOUND"
    message:
        "Finding & fetching significantly associated k-mers for each phenotype..."
    script:
        "../scripts/fetch_significant_kmers.py"

# =========================================================================================================

