# =========================================================================================================
#     Merge reads with k-mers
# =========================================================================================================

## TODO: Find a better way to merge reads here. ##

rule merge_reads:
    input:
        dir= rules.fetch_source_reads.output,
    output:
        merged_r1 = temp("results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_R1.fastq"),
        merged_r2 = temp("results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_R2.fastq")
    shell:
        """
        echo Merging r1 files:
        echo "$(ls {input.dir}/*_reads_with_kmers_R1.fastq)"
        cat {input.dir}/*_reads_with_kmers_R1.fastq > {output.merged_r1}
        echo Merging r2 files:
        echo "$(ls {input.dir}/*_reads_with_kmers_R2.fastq)"
        cat {input.dir}/*_reads_with_kmers_R2.fastq > {output.merged_r2}
        """

# =========================================================================================================
#     Sort merged reads
# =========================================================================================================

rule sort_reads:
    input:
        r1 = "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_R1.fastq",
        r2 = "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_R2.fastq"
    output:
        sorted_r1 = "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_sorted_R1.fastq",
        sorted_r2 = "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_sorted_R2.fastq"
    conda:
        "../envs/align_reads.yaml"
    shell:
        """
        seqkit sort -n {input.r1} > {output.sorted_r1}
        seqkit sort -n {input.r2} > {output.sorted_r2}
        """

# =======================================================================================================
#     Align reads with significant k-mers to the reference genome
# =======================================================================================================

rule align_reads:
    input:
        r1= "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_sorted_R1.fastq",
        r2= "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_sorted_R2.fastq",
        done = "results/fetch_reads_with_kmers/fetch_source_reads.done",
    output:
        out_sam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.sam",
        done = touch("results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.aligning_reads.done")
    params:
        index = "resources/ref/genome/genome",
        extra = config["params"]["bowtie2"]["extra"]
    conda:
        "../envs/align_reads.yaml"
    threads: 
        config["params"]["bowtie2"]["threads"]
    log:
        "logs/align_reads/{phenos_filt}/align_reads_with_kmers.bowtie2.log"
    message:
        "Aligning reads with k-mers to the reference genome..."
    shell:
        """
        bowtie2 -p {threads} --very-sensitive-local {params.extra} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.out_sam} 2> {log}
        """

# =======================================================================================================
#     Filter alignments based on mapping quality
# =======================================================================================================

rule filter_alignment:
    input:
        sam= "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.sam",
    output:
        "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sam",
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

rule align_reads_sam_to_bam:
    input:
        "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sam"
    output:
        temp("results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.bam")
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

rule align_reads_bam_sort:
    input:
        "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.bam"
    output:
        "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam"
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

rule align_reads_bam_index:
    input:
        "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam"
    output:
        "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam.bai"
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
#     Check align_reads 
# =========================================================================================================

def aggregate_input_align_reads(wildcards):
    checkpoint_output = checkpoints.fetch_significant_kmers.get(**wildcards).output[0]
    return expand("results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam.bai",
           phenos_filt=glob_wildcards(os.path.join(checkpoint_output, "{phenos_filt}_kmers_list.txt")).phenos_filt)


rule check_align_reads:
    input:
        aggregate_input_align_reads
    output:
        "results/align_reads_with_kmers/align_reads_with_kmers.done"
    message:
        "Checking if aligning reads with significant k-mers is done..."
    shell:
        """
        touch {output}
        """

# =========================================================================================================