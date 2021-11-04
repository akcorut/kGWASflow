# =================================================================================================
#     Generate kmers_list_paths.txt 
#     # This file contains a list of all individuals/samples k-mers list files and
#     # required to run kmersGWAS.
#     # For more info: https://github.com/voichek/kmersGWAS/blob/master/manual.pdf
# =================================================================================================

rule generate_kmers_list_paths:
    input:
        expand("results/kmers_count/{u.sample_name}/kmers_with_strand", u=samples.itertuples())
    output:
        "results/kmers_list/kmers_list_paths.txt"
    params:
        input_dir=  "results/kmers_count",
        out_dir= lambda wildcards, output: output[0][20]
    shell:
        """
        python scripts/generate_kmers_list_paths.py -i {params.input_dir} -o {params.out_dir}
        """