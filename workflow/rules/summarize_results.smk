# =================================================================================================
#     Generate results table from kmersGWAS results (kmers_pass_threshold_5per)
# =================================================================================================

rule generate_results_table:
    input:
        expand(rules.run_kmers_gwas.output, pheno= pheno_names)
    output:
        res_tab_5per = report(
            "results/kmers_gwas_summary/kmers_gwas_results_table_5per.txt",
            caption="../report/kmers_gwas_results_table.rst",
            category="k-mers GWAS Results Summary",
        ),
        phenos_list_5per= "results/kmers_gwas_summary/phenos_with_sig_kmers_list_5per.txt",
    params:
        in_prefix = "results/kmers_gwas",
        out_prefix = lambda w, output: os.path.dirname(output[0]),
        extra = config["params"]["results_table"]["extra"]
    conda:
        "../envs/kmers_stats.yaml"
    message:
        "Summarizing k-mers GWAS results in a table..."
    shell:
        """
        python scripts/generate_kmers_gwas_results_table.py -dir {params.in_prefix} \
        -out {params.out_prefix} {params.extra}
        """

# =================================================================================================