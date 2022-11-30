# =========================================================================================================
#     Retrieve the significant k-mers (kmers_pass_threshold_5per) from the results table for each phenotype
# =========================================================================================================

checkpoint fetch_significant_kmers:
    input:
        kmers_gwas_res = expand(rules.run_kmers_gwas.output, pheno= pheno_names),
        res_tab = rules.generate_results_table.output.res_tab_5per,
    output:
        directory("results/fetch_kmers/")
    params:
        in_prefix = lambda w, input: os.path.dirname(os.path.dirname(os.path.dirname(input.kmers_gwas_res[0]))), 
        threshold = config["params"]["fetch_kmers"]["threshold"],
        check_file_prefix = lambda w, input: os.path.dirname(input.res_tab)
    log:
        "logs/fetch_kmers/fetch_significant_kmers.log"
    message:
        "Finding & fetching significantly associated k-mers for each phenotype..."
    script:
        "../scripts/fetch_significant_kmers.py"

# =========================================================================================================

