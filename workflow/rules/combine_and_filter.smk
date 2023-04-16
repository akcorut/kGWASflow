# =================================================================================================
#     Generate kmers_list_paths.txt 
#     # This file contains a list of all individuals/samples k-mers list files and
#     # required to run kmersGWAS.
#     # For more info: https://github.com/voichek/kmersGWAS/blob/master/manual.pdf
# =================================================================================================

rule generate_kmers_list_paths:
    input:
        kmers_with_strand = expand("results/kmers_count/{u.sample_name}/kmers_with_strand", u=samples.itertuples()),
        samples_tab = config["samples"]
    output:
        "results/kmers_list/kmers_list_paths.txt"
    params:
        input_dir =  lambda w, input: os.path.dirname(os.path.dirname(input.kmers_with_strand[0])),
        out_dir = lambda wildcards, output: os.path.dirname(output[0])
    log:
        "logs/combine_and_filter/generate_kmers_list_paths.log"
    message:
        "Generating kmers_list_paths.txt..."
    script:
        "../scripts/generate_kmers_list_paths.py"

# =================================================================================================
#     Combine and filter lists of kmers
#     # Source code of 'list_kmers_found_in_multiple_samples': https://github.com/voichek/kmersGWAS
# =================================================================================================

rule combine_and_filter:
    input:
        kmers_count =expand("results/kmers_count/{u.sample_name}/kmers_with_strand", u=samples.itertuples()),
        kmers_list = rules.generate_kmers_list_paths.output,
        kmersGWAS_bin = rules.extract_kmersGWAS.output.kmersGWAS_bin,
    output:
        kmers_to_use = "results/kmers_list/kmers_to_use",
        shareness = "results/kmers_list/kmers_to_use.shareness"
    params:
        kmer_len = config["params"]["kmc"]["kmer_len"],
        mac = config["params"]["kmers_gwas"]["minor_allele_count"],
        min_app = config["params"]["kmers_gwas"]["min_percent_app"]
    conda:
        "../envs/kmers_gwas.yaml"
    log:
        "logs/combine_and_filter/combine_and_filter.log"
    message:
        "Combining the k-mers from each acession/sample into one list and filter the k-mers..."
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib

        {input.kmersGWAS_bin}/list_kmers_found_in_multiple_samples -l {input.kmers_list} -k {params.kmer_len} \
        --mac {params.mac} -p {params.min_app} -o {output} 2> {log}
        """

# =================================================================================================
#     Plot k-mer allele counts
# =================================================================================================

rule plot_kmer_allele_counts:
    input:
        shareness = rules.combine_and_filter.output.shareness
    output:
        kmer_allele_counts_plot = report(
            "results/plots/kmer_allele_counts/kmer_allele_counts.barplot.pdf",
            caption="../report/plot_kmer_allele_counts.rst",
            category="k-mer Counts - Summary Statistics",
            subcategory="k-mer Allele Counts",
        )
    conda:
        "../envs/kmer_stats.yaml"
    message:
        "Plotting the k-mer allele counts..."
    log:
        "logs/plots/kmer_allele_counts/plot_kmer_allele_counts.log"
    script:
        "../scripts/plot_kmer_allele_counts.py"

# =================================================================================================