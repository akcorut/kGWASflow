# =================================================================================================
#     Create the k-mers table
#     # Source code of 'build_kmers_table': https://github.com/voichek/kmersGWAS
# =================================================================================================

rule create_kmers_table:
    input:
        list_paths = rules.generate_kmers_list_paths.output,
        kmers_count =expand("results/kmers_count/{u.sample_name}/kmers_with_strand", u=samples.itertuples()),
        kmers_to_use = rules.combine_and_filter.output.kmers_to_use,
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        "results/kmers_table/kmers_table.table"
    params:
        prefix = lambda w, output: os.path.splitext(output[0])[0],
        kmer_len = config["params"]["kmc"]["kmer_len"],
    conda:
        "../envs/kmers_gwas.yaml"
    log:
        "logs/create_kmers_table/build_kmers_table.log"
    message:
        "Creating the k-mers table that contains the presence absence/pattern of the k-mers..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
        
        {input.kmersGWAS_bin}/build_kmers_table -l {input.list_paths} -k {params.kmer_len} \
        -a {input.kmers_to_use} -o {params.prefix} 2> {log}
        """

# =================================================================================================
#     Convert the k-mers table to PLINK binary file
# =================================================================================================

rule convert_kmers_table_to_plink:
    input:
        kmers_table = "results/kmers_table/kmers_table.table",
        phenotype = lambda wildcards: phenos.loc[(wildcards.pheno), ["pheno_path"]],
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        bed = "results/kmers_table/plink/pheno_{pheno}/kmers_table.presence_absence.0.bed",
        done = touch("results/kmers_table/plink/pheno_{pheno}/convert_kmers_table_to_plink.done"),
    params:
        input_prefix = lambda w, input: os.path.splitext(input.kmers_table)[0],
        output_prefix = lambda w, output: os.path.splitext(output.bed)[0].rsplit(".", 1)[0],
        kmer_len = config["params"]["kmers_table_to_bed"]["kmer_len"],
        max_num_var = config["params"]["kmers_table_to_bed"]["max_num_var"],
        mac = config["params"]["kmers_table_to_bed"]["minor_allele_count"],
        maf = config["params"]["kmers_table_to_bed"]["minor_allele_freq"],
        only_unique = "-u" if config["params"]["kmers_table_to_bed"]["only_unique"] else "",
    conda:
        "../envs/kmers_gwas.yaml"
    log:
        "logs/create_kmers_table/convert_kmers_table_to_plink.pheno_{pheno}.log"
    message:
        "Converting the k-mers table to PLINK binary file..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
        
        {input.kmersGWAS_bin}/kmers_table_to_bed \
        -t {params.input_prefix} \
        -p {input.phenotype} \
        -o {params.output_prefix} \
        -k {params.kmer_len} \
        --maf {params.maf} --mac {params.mac} \
        -b {params.max_num_var} \
        {params.only_unique} \
        >>{log} 2>&1
        """

# =================================================================================================