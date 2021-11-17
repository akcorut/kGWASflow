# =======================================================================================================
#     Align significant k-mers to the reference genome
# =======================================================================================================

rule align_kmers:
    input:
        kmers_list = "results/fetch_kmers/{phenos_filt}_kmers_list.fa",
        index = rules.bowtie2_build.output
    output:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sam"
    params:
        index = "resources/genome",
        m = config["params"]["bowtie"]["max_align_repot"]
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["bowtie"]["threads"]
    log:
        "logs/align_kmers/{phenos_filt}_kmers_align.bowtie2.log"
    message:
        "Aligning signficant k-mers to the reference genome..."
    shell:
        """
        bowtie -p {threads} -a --best --all --strata -m {params.m} \
        -x {params.index} -f {input.kmers_list} --sam {output} 2> {log}
        """

# =======================================================================================================
#     Convert SAM outputs to BAM format
# =======================================================================================================

rule sam_to_bam:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sam"
    output:
        temp("results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.bam")
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    message:
        "Converting SAM files to BAM..."
    shell:
        """
        samtools view -@ {threads} -Sbh {input} > {output}
        """

# =======================================================================================================
#     Sort alignment BAM files
# =======================================================================================================

rule bam_sort:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.bam"
    output:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam"
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    message:
        "Sorting alignment BAM files..."
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        """

# =======================================================================================================
#     Index alignment BAM files
# =======================================================================================================

rule bam_index:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam"
    output:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam.bai"
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    message:
        "Indexing alignment BAM files..."
    shell:
        """
        samtools index -@ {threads} {input}
        """

# =========================================================================================================
#     Check align_kmers
# =========================================================================================================

def aggregate_input_align_kmers(wildcards):
    checkpoint_output = checkpoints.fetch_kmers_from_res_table.get(**wildcards).output[0]
    return expand("results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam.bai",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

rule check_align_kmers:
    input:
        aggregate_input_align_kmers
    output:
        "results/align_kmers/align_kmers.done"
    message:
        "Checking if aligning k-mers is done..."
    shell:
        """
        touch {output}
        """

# =========================================================================================================