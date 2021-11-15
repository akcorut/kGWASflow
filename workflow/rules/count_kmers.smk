# =================================================================================================
#     Create Symlinks
# =================================================================================================

if not config["settings"]["trimming"]["activate"]:
    rule create_symlink:
        input:
            r1 = lambda wildcards: samples.loc[(wildcards.sample, wildcards.library), ["fq1"]],
            r2 = lambda wildcards: samples.loc[(wildcards.sample, wildcards.library), ["fq2"]]
        output:
            r1 = "results/reads/{sample}/{library}_1.fastq",
            r2 = "results/reads/{sample}/{library}_2.fastq",
        message:
            "Creating symbolic links for fastq file..."
        threads: 1
        shell:
            """
            ln -s {input.r1} {output.r1}
            ln -s {input.r2} {output.r2}
            """

# =================================================================================================
#     Generate Lists Of Inputs for KMC
# =================================================================================================

if not config["settings"]["trimming"]["activate"]:
    rule generate_input_lists:
        input:
            r1 = expand("results/reads/{sample}/{library}_1.fastq", zip,
                        sample=sample_names,
                        library=library_names),
            r2 = expand("results/reads/{sample}/{library}_2.fastq", zip,
                        sample=sample_names,
                        library=library_names)
        output:
            "results/reads/{sample}/input_files.txt"
        params:
            prefix = get_input_path_for_generate_input_lists()
        message:
            "Generating input list files..."
        shell:
            "python scripts/generate_input_lists.py -i {params.prefix}"
else:
    rule generate_input_lists:
        output:
            "results/trimmed/{sample}/input_files.txt"
        params:
            prefix = get_input_path_for_generate_input_lists()
        message:
            "Generating input list files..."
        shell:
            "python scripts/generate_input_lists.py -i {params.prefix}"

# =================================================================================================
#     KMC with canonization
# =================================================================================================

rule kmc_canonical:
    input:
        rules.generate_input_lists.output
    output:
        kmc_suf = temp("results/kmers_count/{sample}/output_kmc_canon.kmc_suf"),
        kmc_pre = temp("results/kmers_count/{sample}/output_kmc_canon.kmc_pre")
    params:
        # outdir_prefix = "results/kmers_count/{sample}/",
        outdir_prefix = lambda wildcards, output: output[0][:-24],
        outfile_prefix = lambda wildcards, output: output[0][:-8],
        kmer_len = config["params"]["kmc"]["kmer_len"],
        count_threshold = config["params"]["kmc"]["count_thresh"],
        extra = config["params"]["kmc"]["extra"]
    log:
        "logs/kmc/{sample}/kmc_canon.log"
    conda:
        "../envs/kmc.yaml"
    threads:
        config["params"]["kmc"]["threads"]
    message:
        "Running KMC to count k-mers with canonization..."
    shell:
        """
        kmc -t{threads} {params.extra} -v -k{params.kmer_len} -ci{params.count_threshold} \
        @{input} {params.outfile_prefix} {params.outdir_prefix} \
        1> {params.outdir_prefix}kmc_canon.1 2> {params.outdir_prefix}kmc_canon.2 \
        > {log}
        """

# =================================================================================================
#     KMC without canonization
# =================================================================================================

rule kmc_non_canonical:
    input:
        rules.generate_input_lists.output
    output:
        kmc_suf = temp("results/kmers_count/{sample}/output_kmc_all.kmc_suf"),
        kmc_pre = temp("results/kmers_count/{sample}/output_kmc_all.kmc_pre")
    params:
        # outdir_prefix = "results/kmers_count/{sample}/",
        outdir_prefix = lambda wildcards, output: output[0][:-22],
        outfile_prefix = lambda wildcards, output: output[0][:-8],
        kmer_len = config["params"]["kmc"]["kmer_len"],
        extra = config["params"]["kmc"]["extra"]
    log:
        "logs/kmc/{sample}/kmc_all.log"
    conda:
        "../envs/kmc.yaml"
    threads:
        config["params"]["kmc"]["threads"]
    message:
        "Running KMC to count k-mers without canonization..."
    shell:
        """
        kmc -t{threads} -v -k{params.kmer_len} -ci0 -b \
        @{input} {params.outfile_prefix} {params.outdir_prefix} \
        1> {params.outdir_prefix}kmc_all.1 2> {params.outdir_prefix}kmc_all.2 \
        > {log}
        """

# =================================================================================================
#     Combine the outputs from the two KMC run
#     # Source code of 'kmers_add_strand_information': https://github.com/voichek/kmersGWAS
# =================================================================================================

rule merge_kmers:
    input:
        kmc_canon_suf = rules.kmc_canonical.output.kmc_suf,
        kmc_canon_pre = rules.kmc_canonical.output.kmc_pre,
        kmc_all_suf = rules.kmc_non_canonical.output.kmc_suf,
        kmc_all_pre = rules.kmc_non_canonical.output.kmc_pre,
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        "results/kmers_count/{sample}/kmers_with_strand"
    params:
        prefix_kmc_canon = "results/kmers_count/{sample}/output_kmc_canon",
        prefix_kmc_all = "results/kmers_count/{sample}/output_kmc_all",
        kmer_len = config["params"]["kmc"]["kmer_len"]
    conda:
        "../envs/kmers_gwas.yaml"
    log:
        "logs/kmc/{sample}/add_strand.log.out"
    message:
        "Merging outputs from two KMC k-mers counting results into one list for each sample/individual..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib

        {input.kmersGWAS_bin}/kmers_add_strand_information -c {params.prefix_kmc_canon} -n {params.prefix_kmc_all} -k {params.kmer_len} -o {output} > {log}
        """

# =================================================================================================
#     Plot Useful Stats from k-mers Count
# =================================================================================================

rule kmers_stats:
    input:
        expand("logs/kmc/{u.sample_name}/kmc_all.log", u=samples.itertuples()),
        expand("logs/kmc/{u.sample_name}/kmc_canon.log", u=samples.itertuples()),
    output:
        kmc_canon_joint_plot = report(
            "results/plots/kmers_count/kmc_canon_total_reads_vs_unique_kmers.joint_plot.pdf",
            caption="../report/plot_joint_kmc_canon.rst",
            category="k-mers Count Stats",
        ),
        kmc_all_joint_plot = report(
            "results/plots/kmers_count/kmc_all_total_reads_vs_unique_kmers.joint_plot.pdf",
            caption="../report/plot_joint_kmc_all.rst",
            category="k-mers Count Stats",
        ),
        kmc_canon_plot = report(
            "results/plots/kmers_count/kmc_canon_total_reads_vs_unique_kmers.scatter_plot.pdf",
            caption="../report/plot_scatter_kmc_canon.rst",
            category="k-mers Count Stats",
        ),
        kmc_all_plot = report(
            "results/plots/kmers_count/kmc_all_total_reads_vs_unique_kmers.scatter_plot.pdf",
            caption="../report/plot_scatter_kmc_all.rst",
            category="k-mers Count Stats",
        ),
        kmc_canon_stats = report(
            "results/tables/kmers_count/kmc_canon.stats.tsv",
            caption="../report/kmc_canon_count_stats.rst",
            category="k-mers Count Stats",
        ),
         kmc_all_stats = report(
            "results/tables/kmers_count/kmc_all.stats.tsv",
            caption="../report/kmc_all_count_stats.rst",
            category="k-mers Count Stats",
        )
    params:
        input_path= "logs/kmc",
        out_table_path= lambda w, output: os.path.dirname(output.kmc_all_stats),
        out_plot_path= lambda w, output: os.path.dirname(output.kmc_all_plot)
    conda:
        "../envs/kmers_stats.yaml"
    message:
        "Gathering useful stats from KMC runs and ploting the results..."
    shell:
        """
        python scripts/plot_kmers_stats.py -i {params.input_path} -o1 {params.out_table_path} -o2 {params.out_plot_path}
        """

# =================================================================================================