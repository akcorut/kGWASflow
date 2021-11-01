# =================================================================================================
#     Create Symlinks
# =================================================================================================

rule create_symlink:
    input:
        r1= lambda wildcards: samples.loc[(wildcards.sample, wildcards.library), ["fq1"]],
        r2= lambda wildcards: samples.loc[(wildcards.sample, wildcards.library), ["fq1"]]
    output:
        r1 =temp("results/reads/{sample}/{library}_1.fastq.gz"),
        r2 =temp("results/reads/{sample}/{library}_2.fastq.gz"),
    message:
        "Create a symbolic link"
    threads: 1
    wildcard_constraints:
        read= [1,2]
    shell:
        "ln -s {input} {output}"

# =================================================================================================
#     Generate Lists Of Inputs for KMC
# =================================================================================================

rule generate_input_lists:
    input:
        get_input_path_for_input_list
    output:
        
    shell:
        "python scripts/generate_input_lists.py -i {input}"