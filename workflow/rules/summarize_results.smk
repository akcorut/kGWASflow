# =================================================================================================
#     Generate results table from kmersGWAS results
# =================================================================================================

rule generate_results_summary:
    input:
        kmers_gwas_res = "results/kmers_gwas/{pheno}/kmers/pass_threshold_5per"
    output:
        res_sum_5per = report(
            "results/tables/kmers_gwas_summary/{pheno}/kmers_gwas.results_summary.5per.txt",
            caption="../report/kmers_gwas_results_summary.rst",
            category= "k-mers GWAS Results - Summary Statistics",
            subcategory="{pheno}",
        ),
        res_sum_10per = report(
            "results/tables/kmers_gwas_summary/{pheno}/kmers_gwas.results_summary.10per.txt",
            caption="../report/kmers_gwas_results_summary.rst",
            category= "k-mers GWAS Results - Summary Statistics",
            subcategory="{pheno}",
        ),
        pval_hist_dir = report(
            directory("results/plots/pval_hist/{pheno}"),
            patterns= ["{pheno}.kmers_pass_first_filter.pval_hist.pdf"],
            caption= "../report/kmers_gwas_results_summary.rst",
            category= "k-mers GWAS Results - Summary Statistics",
            subcategory= "{pheno}"),
        done = touch("results/tables/kmers_gwas_summary/{pheno}/generate_results_summary.done")
    params:
        pheno = "{pheno}",
        kmers_number = config["params"]["kmers_gwas"]["kmers_number"],
    conda:
        "../envs/kmer_stats.yaml"
    log:
        "logs/tables/kmers_gwas_summary/{pheno}/generate_kmers_gwas_results_summary.log"
    message:
        "Summarizing k-mers GWAS results in a table..."
    script:
        "../scripts/generate_kmers_gwas_summary.py"

# =================================================================================================