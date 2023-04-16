# =================================================================================================
#     Create Symlinks
# =================================================================================================

## TODO: Find a better solution here for the optional output (.gz or not .gz) ##

R1_OUT=""
R2_OUT=""
if ends_with_gz(samples):
    R1_OUT= "results/reads/{sample}/{library}_1.fastq.gz"
    R2_OUT= "results/reads/{sample}/{library}_2.fastq.gz"
else:
    R1_OUT= "results/reads/{sample}/{library}_1.fastq"
    R2_OUT= "results/reads/{sample}/{library}_2.fastq"


if not config["settings"]["trimming"]["activate"]:
    rule create_symlink:
        input:
            fastqs= get_fastqs
        output:
            r1 = R1_OUT,
            r2 = R2_OUT
        log:
            log1="logs/count_kmers/create_symlink/{sample}/{sample}_{library}_1.fastq.create_symlink.log",
            log2="logs/count_kmers/create_symlink/{sample}/{sample}_{library}_2.fastq.create_symlink.log",
        message:
            "Creating symbolic links for fastq files..."
        threads: 1
        shell:
            """
            echo Working on fastq files: {input.fastqs}
            echo Symlink -fastq1: {input.fastqs[0]} to {output.r1}
            ln -rs {input.fastqs[0]} {output.r1} 2> {log.log1}
            echo Symlink -fastq2: {input.fastqs[1]} to {output.r2}
            ln -rs {input.fastqs[1]} {output.r2} 2> {log.log2}
            """

# =================================================================================================
#     Generate Lists Of Inputs for KMC
# =================================================================================================

if not config["settings"]["trimming"]["activate"]:
    rule generate_input_lists:
        input:
            r1 = expand(R1_OUT, zip,
                        sample=sample_names,
                        library=library_names),
            r2 = expand(R2_OUT, zip,
                        sample=sample_names,
                        library=library_names),
            qc= "results/qc/multiqc.html"
        output:
            "results/reads/{sample}/input_files.txt"
        params:
            prefix = lambda w, output: os.path.dirname(os.path.dirname(output[0]))
        priority: 60
        log:
            "logs/count_kmers/generate_input_lists/{sample}/{sample}_generate_input_lists.log"
        message:
            "Generating input list files..."
        script:
            "../scripts/generate_input_lists.py"

if config["settings"]["trimming"]["activate"]:
    rule generate_input_lists:
        input:
            r1 = expand(rules.cutadapt_pe.output.fastq1, zip,
                        sample=sample_names,
                        library=library_names),
            r2 = expand(rules.cutadapt_pe.output.fastq2, zip,
                        sample=sample_names,
                        library=library_names),
            qc= "results/qc/multiqc.html"
        output:
            "results/trimmed/{sample}/input_files.txt"
        params:
            prefix = lambda w, output: os.path.dirname(os.path.dirname(output[0]))
        priority: 60
        log:
            "logs/count_kmers/generate_input_lists/{sample}/{sample}_generate_input_lists.log"
        message:
            "Generating input list files..."
        script:
            "../scripts/generate_input_lists.py"

# =================================================================================================
#     KMC with canonization
# =================================================================================================

rule kmc_canonical:
    input:
        rules.generate_input_lists.output
    output:
        kmc_suf = temp("results/kmers_count/{sample}/output_kmc_canon.kmc_suf"),
        kmc_pre = temp("results/kmers_count/{sample}/output_kmc_canon.kmc_pre"),
        done= touch("results/kmers_count/{sample}/kmc_canonical.done")
    params:
        outdir_prefix = lambda w, output: os.path.dirname(output.kmc_suf),
        outfile_prefix = lambda w, output: os.path.splitext(output.kmc_suf)[0],
        kmer_len = config["params"]["kmc"]["kmer_len"],
        count_threshold = config["params"]["kmc"]["count_thresh"],
        extra = config["params"]["kmc"]["extra"]
    log:
        "logs/count_kmers/kmc/{sample}/kmc_canon.log"
    conda:
        "../envs/kmc.yaml"
    threads:
        config["params"]["kmc"]["threads"]
    priority: 50
    retries: 3
    message:
        "Running KMC to count k-mers with canonization..."
    shell:
        """
        kmc -t{threads} {params.extra} -v -k{params.kmer_len} -ci{params.count_threshold} \
        @{input} {params.outfile_prefix} {params.outdir_prefix} \
        1> {params.outdir_prefix}/kmc_canon.1 2> {params.outdir_prefix}/kmc_canon.2 \
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
        kmc_pre = temp("results/kmers_count/{sample}/output_kmc_all.kmc_pre"),
        done= touch("results/kmers_count/{sample}/kmc_non-canonical.done")
    params:
        outdir_prefix = lambda w, output: os.path.dirname(output.kmc_suf),
        outfile_prefix = lambda w, output: os.path.splitext(output.kmc_suf)[0],
        kmer_len = config["params"]["kmc"]["kmer_len"],
        extra = config["params"]["kmc"]["extra"]
    log:
        "logs/count_kmers/kmc/{sample}/kmc_all.log"
    conda:
        "../envs/kmc.yaml"
    threads:
        config["params"]["kmc"]["threads"]
    priority: 40
    retries: 3
    message:
        "Running KMC to count k-mers without canonization..."
    shell:
        """
        kmc -t{threads} -v -k{params.kmer_len} -ci0 -b \
        @{input} {params.outfile_prefix} {params.outdir_prefix} \
        1> {params.outdir_prefix}/kmc_all.1 2> {params.outdir_prefix}/kmc_all.2 \
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
        kmc_canon_done = rules.kmc_canonical.output.done,
        kmc_all_done = rules.kmc_non_canonical.output.done,
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        strand= temp("results/kmers_count/{sample}/kmers_with_strand"),
        done= touch("results/kmers_count/{sample}/kmers_add_strand_information.done")
    params:
        prefix_kmc_canon = lambda w, input: os.path.splitext(input.kmc_canon_suf)[0],
        prefix_kmc_all = lambda w, input: os.path.splitext(input.kmc_all_suf)[0],
        kmer_len = config["params"]["kmc"]["kmer_len"]
    conda:
        "../envs/kmers_gwas.yaml"
    threads:
        config["params"]["merge_kmers"]["threads"]
    log:
        "logs/count_kmers/kmc/{sample}/add_strand.log.out"
    retries: 3
    message:
        "Merging outputs from two KMC k-mers counting results into one list for each sample/individual..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib

        {input.kmersGWAS_bin}/kmers_add_strand_information -c {params.prefix_kmc_canon} -n {params.prefix_kmc_all} -k {params.kmer_len} -o {output.strand} > {log}
        """

# =================================================================================================
#     Plot Useful Stats from k-mers Count
# =================================================================================================

rule kmer_stats:
    input:
        kmc_all = expand("logs/count_kmers/kmc/{u.sample_name}/kmc_all.log", u=samples.itertuples()),
        kmc_canon = expand("logs/count_kmers/kmc/{u.sample_name}/kmc_canon.log", u=samples.itertuples()),
    output:
        kmc_canon_joint_plot = report(
            "results/plots/kmer_stats/kmc_canon.total_reads_vs_unique_kmers.joint_plot.pdf",
            caption="../report/plot_joint_kmc_canon.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="Total Reads vs. Unique k-mers",
        ),
        kmc_all_joint_plot = report(
            "results/plots/kmer_stats/kmc_all.total_reads_vs_unique_kmers.joint_plot.pdf",
            caption="../report/plot_joint_kmc_all.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="Total Reads vs. Unique k-mers",
        ),
        kmc_canon_plot = report(
            "results/plots/kmer_stats/kmc_canon.total_reads_vs_unique_kmers.scatter_plot.pdf",
            caption="../report/plot_scatter_kmc_canon.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="Total Reads vs. Unique k-mers",
        ),
        kmc_all_plot = report(
            "results/plots/kmer_stats/kmc_all.total_reads_vs_unique_kmers.scatter_plot.pdf",
            caption="../report/plot_scatter_kmc_all.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="Total Reads vs. Unique k-mers",
        ),
        kmc_all_kmer_dist_plot = report(
            "results/plots/kmer_dist/kmc_all.kmer_dist_hist.pdf",
            caption="../report/plot_kmer_dist_hist.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="k-mer Counts Distribution",
        ),
        kmc_canon_kmer_dist_plot = report(
            "results/plots/kmer_dist/kmc_canon.kmer_dist_hist.pdf",
            caption="../report/plot_kmer_dist_hist.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="k-mer Counts Distribution",
        ),
        kmc_canon_stats = report(
            "results/tables/kmer_stats/kmc_canon.kmer_stats.tsv",
            caption="../report/kmc_canon_count_stats.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="Total Reads vs. Unique k-mers",
        ),
         kmc_all_stats = report(
            "results/tables/kmer_stats/kmc_all.kmer_stats.tsv",
            caption="../report/kmc_all_count_stats.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="Total Reads vs. Unique k-mers",
        )
    params:
        input_path= lambda w, input: os.path.dirname(os.path.dirname(input.kmc_canon[0]))
    log:
        "logs/plots/kmer_stats/plot_kmer_stats.log"
    conda:
        "../envs/kmer_stats.yaml"
    message:
        "Gathering useful stats from KMC runs and ploting the results..."
    script:
        "../scripts/plot_kmer_stats.py"

# =================================================================================================