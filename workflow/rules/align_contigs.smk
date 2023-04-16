# =========================================================================================================
#     Align contigs to reference genome - Minimap2
# =========================================================================================================

rule align_contigs_minimap:
    input:
        contigs= "results/assemble_reads_with_kmers/{phenos_filt}/assembly/contigs.fasta",
        ref_gen= "resources/ref/genome/genome.fasta",
    output:
        "results/align_contigs/{phenos_filt}/alignment/contigs_aligned.sam"
    conda:
        "../envs/align_contigs.yaml"
    threads: 
        config["params"]["minimap2"]["threads"]
    log:
        "logs/align_contigs/{phenos_filt}/align_contigs.minimap2.log"
    message:
        "Aligning assembled contigs of reads with associated k-mers to the reference genome..."
    shell:
        """
        minimap2 -t {threads} -a {input.ref_gen} {input.contigs} > {output} 2> {log}
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
    log:
        "logs/align_contigs/{phenos_filt}/align_contigs.filter_alignment.log"
    message:
        "Filtering alignment results based on mapping quality..."
    shell:
        """
        awk -v s={params.min_mapping_score} '$5 > s || $1 ~ /^@/' {input} > {output} 2> {log}
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
        "../envs/align_contigs.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_contigs/{phenos_filt}/align_contigs.sam_to_bam.log"
    message:
        "Converting SAM files to BAM..."
    shell:
        """
        samtools view -@ {threads} -Sbh {input} > {output} 2> {log}
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
        "../envs/align_contigs.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_contigs/{phenos_filt}/align_contigs.bam_sort.log"
    message:
        "Sorting alignment BAM files..."
    shell:
        """
        samtools sort -@ {threads} {input} -o {output} 2> {log}
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
        "../envs/align_contigs.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_contigs/{phenos_filt}/align_contigs.bam_index.log"
    message:
        "Indexing alignment BAM files..."
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """

# =======================================================================================================
#     Convert alignment BAM files to BED format
# =======================================================================================================

rule align_contigs_bam_to_bed:
    input:
        bam="results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam",
        bai="results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam.bai"
    output:
        "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bed"
    conda:
        "../envs/bedtools.yaml"
    threads:
        config["params"]["bedtools"]["threads"]
    log:
        "logs/align_contigs/{phenos_filt}/align_contigs.bam_to_bed.log"
    message:
        "Converting alignment BAM files to BED..."
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output} 2> {log}
        """

# =========================================================================================================
#     Aggregate align_contigs outputs
# =========================================================================================================

rule aggregate_align_contigs:
    input:
        aggregate_input_align_contigs
    output:
        "results/align_contigs/align_contigs.done"
    log:
        "logs/align_contigs/aggregate_align_contigs.log"
    message:
        "Checking if aligning contigs is done..."
    shell:
        """
        touch {output} 2> {log}
        """

# =========================================================================================================