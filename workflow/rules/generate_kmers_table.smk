# =================================================================================================
#     Create the k-mers table
# =================================================================================================

rule create_kmers_table:
    input:
        list_paths = rules.generate_kmers_list_paths.output,
        kmers_to_use = rules.combine_and_filter.output.kmers_to_use,
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        "results/kmers_table/kmers_table.table"
    params:
        prefix = lambda wildcards, output: output[0][:-6],
        kmer_len = config["params"]["kmc"]["kmer_len"],
    conda:
        "../envs/kmers_gwas.yaml"
    message:
        "Creating the k-mers table that contains the presence absence/pattern of the k-mers..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
        
        {input.kmersGWAS_bin}/build_kmers_table -l {input.list_paths} -k {params.kmer_len} -a {input.kmers_to_use} -o {params.prefix}
        """

# =================================================================================================
#     TODO Convert the k-mers table to PLINK binary file
# =================================================================================================