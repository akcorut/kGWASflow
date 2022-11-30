# =================================================================================================
#     Generate results table from kmersGWAS results (kmers_pass_threshold_5per)
# =================================================================================================

rule generate_results_table:
    input:
        kmers_gwas_res = expand(rules.run_kmers_gwas.output, pheno= pheno_names)
    output:
        res_tab_5per = report(
            "results/kmers_gwas_summary/kmers_gwas_results_table_5per.txt",
            caption="../report/kmers_gwas_results_table.rst",
            category="k-mers GWAS Results Summary",
        ),
        phenos_list_5per= "results/kmers_gwas_summary/phenos_with_sig_kmers_list_5per.txt",
    params:
        in_prefix = lambda w, input: os.path.dirname(os.path.dirname(os.path.dirname(input.kmers_gwas_res[0]))),
        out_prefix = lambda w, output: os.path.dirname(output[0]),
    conda:
        "../envs/kmers_stats.yaml"
    log:
        "logs/kmers_gwas_summary/generate_kmers_gwas_results_table.log"
    message:
        "Summarizing k-mers GWAS results in a table..."
    script:
        "../scripts/generate_kmers_gwas_results_table.py"

# =================================================================================================