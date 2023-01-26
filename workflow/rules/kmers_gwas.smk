# =================================================================================================
#     Run kmers-based GWAS
#     # This method adapted from Voichek et al. (2020)
#     # Reference: https://www.nature.com/articles/s41588-020-0612-7
# =================================================================================================

rule run_kmers_gwas:
    input:
        kmers_gwas_py = rules.extract_kmersGWAS.output.kmersGWAS_py,
        phenotype = lambda wildcards: phenos.loc[(wildcards.pheno), ["pheno_path"]],
        kmers_table = rules.create_kmers_table.output,
        kinship_matrix = "results/kmers_table/kmers_table.kinship",
    output:
        "results/kmers_gwas/{pheno}/kmers/pass_threshold_5per"
    params:
        kmers_tab_prefix = lambda w, input: os.path.splitext(input.kinship_matrix)[0],
        out_prefix = lambda wildcards, output: output[0][:-26],
        kmer_len = config["params"]["kmc"]["kmer_len"],
        min_data_points = config["params"]["kmers_gwas"]["min_data_points"],
        mac = config["params"]["kmers_gwas"]["minor_allele_count"],
        maf = config["params"]["kmers_gwas"]["minor_allele_freq"],
        kmers_number = config["params"]["kmers_gwas"]["kmers_number"],
        n_permutations = config["params"]["kmers_gwas"]["n_permutations"],
        extra = config["params"]["kmers_gwas"]["extra"]
    threads:
        config["params"]["kmers_gwas"]["threads"]
    log:
        "logs/kmers_gwas/{pheno}/{pheno}_kmers_gwas_run.log"
    conda:
        "../envs/kmers_gwas_py2.yaml"
    message:
        "Running k-mers based GWAS on {input.phenotype}..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
        
        rm -r {params.out_prefix}

        python2 {input.kmers_gwas_py} {params.extra} \
        --min_data_points {params.min_data_points} \
        --pheno {input.phenotype} \
        --kmers_table {params.kmers_tab_prefix} \
        --kmers_number {params.kmers_number} \
        --permutations {params.n_permutations} \
        --maf {params.maf} --mac {params.mac} \
        -l {params.kmer_len} -p {threads} \
        --outdir {params.out_prefix} >>{log} 2>&1
        """

# =================================================================================================