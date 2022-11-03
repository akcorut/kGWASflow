# =======================================================================================================
#     Align significant k-mers to the reference genome
# =======================================================================================================

if config["settings"]["align_kmers"]["use_bowtie"]:
    rule align_kmers:
        input:
            kmers_list = "results/fetch_kmers/{phenos_filt}_kmers_list.fa",
            index = rules.bowtie_build.output
        output:
            "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sam"
        params:
            index = "resources/ref/genome/bowtie_index/genome",
            extra = config["params"]["bowtie"]["extra"]
        conda:
            "../envs/align_kmers.yaml"
        threads:
            config["params"]["bowtie"]["threads"]
        log:
            "logs/align_kmers/{phenos_filt}_kmers_align.bowtie.log"
        message:
            "Aligning signficant k-mers to the reference genome..."
        shell:
            """
            bowtie -p {threads} -a --best --strata {params.extra} \
            -x {params.index} -f {input.kmers_list} --sam {output} 2> {log}
            """


if config["settings"]["align_kmers"]["use_bowtie2"]:
    rule align_kmers:
        input:
            kmers_list = "results/fetch_kmers/{phenos_filt}_kmers_list.fa",
            index = rules.bowtie2_build.output
        output:
            "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sam"
        params:
            index = "resources/ref/genome/bowtie2_index/genome",
            extra = config["params"]["bowtie2"]["extra"]
        conda:
            "../envs/align_kmers.yaml"
        threads:
            config["params"]["bowtie2"]["threads"]
        log:
            "logs/align_kmers/{phenos_filt}_kmers_align.bowtie2.log"
        message:
            "Aligning signficant k-mers to the reference genome..."
        shell:
            """
            bowtie2 -p {threads} {params.extra} \
            -x {params.index} -f {input.kmers_list} -S {output} 2> {log}
            """

# =======================================================================================================
#     Convert SAM outputs to BAM format
# =======================================================================================================

rule align_kmers_sam_to_bam:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sam"
    output:
        temp("results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.bam")
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_kmers/{phenos_filt}_kmers_align.sam_to_bam.log"
    message:
        "Converting SAM files to BAM..."
    shell:
        """
        samtools view -@ {threads} -Sbh {input} > {output} 2> {log}
        """

# =======================================================================================================
#     Sort alignment BAM files
# =======================================================================================================

rule align_kmers_bam_sort:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.bam"
    output:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam"
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_kmers/{phenos_filt}_kmers_align.bam_sort.log"
    message:
        "Sorting alignment BAM files..."
    shell:
        """
        samtools sort -@ {threads} {input} -o {output} 2> {log}
        """

# =======================================================================================================
#     Index alignment BAM files
# =======================================================================================================

rule align_kmers_bam_index:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam"
    output:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sorted.bam.bai"
    conda:
        "../envs/align_kmers.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_kmers/{phenos_filt}_kmers_align.bam_index.log"
    message:
        "Indexing alignment BAM files..."
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """

# =================================================================================================
#     Generate Manhattan plots
# =================================================================================================

rule plot_manhattan:
    input:
        "results/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.sam"
    output:
        manhattan_plot = report(
            "results/plots/manhattan/align_kmers/{phenos_filt}/{phenos_filt}_kmers_alignment.manhattan_plot.pdf",
            caption="../report/plot_manhattan.rst",
            category="k-mers GWAS Results",
        )
    conda:
        "../envs/plot_manhattan.yaml"
    threads:
        config["params"]["plot_manhattan"]["threads"]
    log:
        "logs/plots/manhattan/align_kmers/{phenos_filt}_kmers_alignment.plot_manhattan.log"
    message:
        "Generating manhattan plot from k-mers alignment results..."
    script:
        "../scripts/plot_manhattan.R"

# =========================================================================================================
#     Aggregate align_kmers outputs
# =========================================================================================================

rule aggregate_align_kmers:
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