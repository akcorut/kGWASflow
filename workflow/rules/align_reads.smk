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
    log:
        "logs/align_reads/{phenos_filt}/sort_reads.seqkit.log"
    message:
        "Sorting reads..."
    shell:
        """
        seqkit sort -n {input.r1} > {output.sorted_r1} 2> {log}
        seqkit sort -n {input.r2} > {output.sorted_r2} 2> {log}
        """

# =======================================================================================================
#     Align reads with significant k-mers to the reference genome
# =======================================================================================================

rule align_reads:
    input:
        r1= "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_sorted_R1.fastq",
        r2= "results/fetch_reads_with_kmers/{phenos_filt}/reads_with_kmers_from_all_acc_sorted_R2.fastq",
        index = rules.bowtie2_build.output,
        fetch_source_reads = "results/fetch_reads_with_kmers/fetch_source_reads.done",
    output:
        out_sam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.sam",
    params:
        index = lambda w, input: os.path.splitext(input.index[0])[0].split('.')[0],
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
    log:
        "logs/align_reads/{phenos_filt}/filter_alignment.log"
    message:
        "Filtering alignment results based on mapping quality..."
    shell:
        """
        awk -v s={params.min_mapping_score} '$5 > s || $1 ~ /^@/' {input} > {output} 2> {log}
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
    log:
        "logs/align_reads/{phenos_filt}/align_reads.sam_to_bam.log"
    message:
        "Converting SAM files to BAM..."
    shell:
        """
        samtools view -@ {threads} -Sbh {input} > {output} 2> {log}
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
    log:
        "logs/align_reads/{phenos_filt}/align_reads.bam_sort.log"
    message:
        "Sorting alignment BAM files..."
    shell:
        """
        samtools sort -@ {threads} {input} -o {output} 2> {log}
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
    log:
        "logs/align_reads/{phenos_filt}/align_reads.bam_index.log"
    message:
        "Indexing alignment BAM files..."
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """

# =========================================================================================================
#     Aggregate align_reads outputs
# =========================================================================================================

rule aggregate_align_reads:
    input:
        aggregate_input_align_reads
    output:
        "results/align_reads_with_kmers/align_reads_with_kmers.done"
    log:
        "logs/align_reads/aggregate_align_reads.log"
    message:
        "Checking if aligning reads with significant k-mers is done..."
    shell:
        """
        touch {output} 2> {log}
        """

# =========================================================================================================