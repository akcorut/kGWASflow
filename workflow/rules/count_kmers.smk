# =================================================================================================
#     Create Symlinks
# =================================================================================================

if not config["settings"]["trimming"]["activate"]:
    rule create_symlink:
        input:
            r1= lambda wildcards: samples.loc[(wildcards.sample, wildcards.library), ["fq1"]],
            r2= lambda wildcards: samples.loc[(wildcards.sample, wildcards.library), ["fq2"]]
        output:
            r1= "results/reads/{sample}/{library}_1.fastq",
            r2= "results/reads/{sample}/{library}_2.fastq",
        message:
            "Create a symbolic link"
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
        shell:
            "python scripts/generate_input_lists.py -i {params.prefix}"
else:
    rule generate_input_lists:
        output:
            "results/trimmed/{sample}/input_files.txt"
        params:
            prefix = get_input_path_for_generate_input_lists()
        shell:
            "python scripts/generate_input_lists.py -i {params.prefix}"

# =================================================================================================
#     KMC with canonization
# =================================================================================================

