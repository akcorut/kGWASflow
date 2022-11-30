# =================================================================================================
#     Generate kinship matrix (k-mers based)
#     # Source code of 'emma_kinship_kmers': https://github.com/voichek/kmersGWAS
# =================================================================================================

if config["settings"]["kmers_gwas"]["use_kmers_kinship"]["activate"]:
    rule generate_kmers_kinship_matrix:
        input:
            kmers_table= rules.create_kmers_table.output,
            kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
        output:
            "results/kmers_table/kmers_table.kinship"
        params:
            prefix = lambda wildcards, output: output[0][:-8],
            kmer_len = config["params"]["kmc"]["kmer_len"],
            mac = config["params"]["kmers_gwas"]["minor_allele_count"],
            maf = config["params"]["kmers_gwas"]["minor_allele_freq"],
        conda:
            "../envs/kmers_gwas.yaml"
        log:
            "logs/generate_kmers_kinship_matrix/emma_kinship_kmers.log"
        message:
            "Calculating the kinship matrix from k-mers..."
        shell:
            """
            export LD_LIBRARY_PATH=$CONDA_PREFIX/lib

            {input.kmersGWAS_bin}/emma_kinship_kmers -t {params.prefix} -k {params.kmer_len} \
            --maf {params.maf} > {output} 2> {log}
            """

# =================================================================================================
#     Generate kinship matrix (SNPs based)
#     # Source code of 'emma_kinship': https://github.com/voichek/kmersGWAS
# =================================================================================================

if not config["settings"]["kmers_gwas"]["use_kmers_kinship"]["activate"]:
    rule generate_snps_kinship_matrix:
        input:
            snps_plink = config["settings"]["kmers_gwas"]["use_snps_kinship"]["snps_plink"],
            kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin
        output:
            "results/kmers_table/kmers_table.kinship"
        params:
            prefix = get_plink_prefix(),
        conda:
            "../envs/kmers_gwas.yaml"
        log:
            "logs/generate_kmers_kinship_matrix/emma_kinship.log"
        message:
            "Generating the kinship matrix from SNPs using {input.snps_plink}..."
        shell:
            """
            export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
            
            {input.kmersGWAS_bin}/emma_kinship {params.prefix} > {output} 2> {log}
            """

# =================================================================================================