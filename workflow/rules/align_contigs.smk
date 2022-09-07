# =========================================================================================================
#     Align contigs to reference genome - Minimap2
# =========================================================================================================

rule align_contigs_minimap:
    input:
        contigs= "results/assemble_reads_with_kmers/{phenos_filt}/assembly/contigs.fasta",
        ref_gen= "resources/genome.fasta",
    output:
        "results/align_contigs/{phenos_filt}/alignment/contigs_aligned.sam"
    conda:
        "../envs/align_contigs.yaml"
    threads: 
        config["params"]["minimap2"]["threads"]
    shell:
        """
        minimap2 -t {threads} -a {input.ref_gen} {input.contigs} > {output}
        """

# =======================================================================================================
#     Filter alignments based on mapping quality
# =======================================================================================================

rule filter_alignment_contigs:
    input:
        sam= "results/align_contigs/{phenos_filt}/alignment/contigs_aligned.sam",
    output:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sam",
    params:
        min_mapping_score= config["params"]["filter_alignment"]["min_map_score"]
    message:
        "Filtering alignment results based on mapping quality..."
    shell:
        """
        awk -v s={params.min_mapping_score} '$5 > s || $1 ~ /^@/' {input} > {output}
        """

# =======================================================================================================
#     Convert SAM outputs to BAM format
# =======================================================================================================

rule align_contigs_sam_to_bam:
    input:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sam"
    output:
        temp("results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.bam")
    conda:
        "../envs/align_reads.yaml"
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

rule align_contigs_bam_sort:
    input:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.bam"
    output:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam"
    conda:
        "../envs/align_reads.yaml"
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

rule align_contigs_bam_index:
    input:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam"
    output:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam.bai"
    conda:
        "../envs/align_reads.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    message:
        "Indexing alignment BAM files..."
    shell:
        """
        samtools index -@ {threads} {input}
        """

# =========================================================================================================
#     Check extract_paired_reads 
# =========================================================================================================

def aggregate_input_align_contigs(wildcards):
    checkpoint_output = checkpoints.fetch_kmers_from_res_table.get(**wildcards).output[0]
    return expand("results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam.bai",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)

rule check_align_contigs:
    input:
        aggregate_input_align_contigs
    output:
        "results/align_contigs/align_contigs.done"
    shell:
        """
        touch {output}
        """