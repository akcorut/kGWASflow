# =========================================================================================================
#     Retrieve the significant k-mers (kmers_pass_threshold_5per) from the results table for each phenotype
# =========================================================================================================

checkpoint fetch_significant_kmers:
    input:
        kmers_gwas_res = expand("results/kmers_gwas/{pheno}/kmers/pass_threshold_5per", pheno = pheno_names),
        res_sum = expand("results/tables/kmers_gwas_summary/{pheno}/generate_results_summary.done", pheno = pheno_names),
    output:
        dir = directory("results/fetch_kmers"),
    params:
        in_prefix = lambda w, input: os.path.dirname(os.path.dirname(input.kmers_gwas_res[0])), 
        threshold = config["params"]["fetch_kmers"]["threshold"],
    log:
        "logs/fetch_kmers/fetch_significant_kmers.log"
    message:
        "Finding & fetching significantly associated k-mers for each phenotype..."
    script:
        "../scripts/fetch_significant_kmers.py"

# =========================================================================================================

