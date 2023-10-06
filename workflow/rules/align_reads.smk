# =========================================================================================================
#     Merge reads with k-mers
# =========================================================================================================

## TODO: Find a better way to merge reads here. ##

rule merge_reads:
    input:
        dir= "results/fetch_reads_with_kmers/{phenos_filt}/individual_reads",
    output:
        merged_r1 = temp("results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.R1.fastq"),
        merged_r2 = temp("results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.R2.fastq"),
        done= touch("results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}.merging_source_reads.done")
    log:
        "logs/align_reads/{phenos_filt}/merge_reads.log"
    shell:
        """
        echo Merging r1 files:
        echo "$(ls {input.dir}/*_reads_with_kmers_R1.fastq)"
        cat {input.dir}/*_reads_with_kmers_R1.fastq > {output.merged_r1} 2> {log}
        echo Merging r2 files:
        echo "$(ls {input.dir}/*_reads_with_kmers_R2.fastq)"
        cat {input.dir}/*_reads_with_kmers_R2.fastq > {output.merged_r2} 2> {log}
        """

# =========================================================================================================
#     Sort merged reads
# =========================================================================================================

rule sort_reads:
    input:
        r1 = "results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.R1.fastq",
        r2 = "results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.R2.fastq"
    output:
        sorted_r1 = "results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.sorted.R1.fastq",
        sorted_r2 = "results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.sorted.R2.fastq",
        done = touch("results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}.sorting_merged_source_reads.done")
    params:
        rename_duplicates = config["params"]["sort_reads"]["rename_dups"]
    conda:
        "../envs/align_reads.yaml"
    log:
        "logs/align_reads/{phenos_filt}/sort_reads.seqkit.log"
    message:
        "Sorting reads..."
    shell:
        """
        if [ {params.rename_duplicates} = True ]; then
            seqkit rename {input.r1} | seqkit sort -n > {output.sorted_r1} 2> {log}
            seqkit rename {input.r2} | seqkit sort -n > {output.sorted_r2} 2>> {log}
        else
            seqkit sort -n {input.r1} > {output.sorted_r1} 2> {log}
            seqkit sort -n {input.r2} > {output.sorted_r2} 2>> {log}
        fi
        """

# =======================================================================================================
#     Align reads with significant k-mers to the reference genome
# =======================================================================================================

rule align_reads:
    input:
        r1= "results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.sorted.R1.fastq",
        r2= "results/fetch_reads_with_kmers/{phenos_filt}/merged_reads/{phenos_filt}_reads_with_kmers.merged.sorted.R2.fastq",
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
        sam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.sam",
    output:
        filtered_sam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sam",
    params:
        min_mapping_score= config["params"]["filter_alignment"]["min_map_score"]
    log:
        "logs/align_reads/{phenos_filt}/filter_alignment.log"
    message:
        "Filtering alignment results based on mapping quality..."
    shell:
        """
        awk -v s={params.min_mapping_score} '$5 > s || $1 ~ /^@/' {input.sam} > {output.filtered_sam} 2> {log}
        """

# =======================================================================================================
#     Convert SAM outputs to BAM format
# =======================================================================================================

rule align_reads_sam_to_bam:
    input:
        sam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sam"
    output:
        bam = temp("results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.bam")
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
        samtools view -@ {threads} -Sbh {input.sam} > {output.bam} 2> {log}
        """

# =======================================================================================================
#     Sort alignment BAM files
# =======================================================================================================

rule align_reads_bam_sort:
    input:
        bam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.bam"
    output:
        sorted_bam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam",
        bam_sort_done = touch("results/align_reads_with_kmers/{phenos_filt}/bam_sort.done")
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
        samtools sort -@ {threads} {input.bam} -o {output.sorted_bam} 2> {log}
        """

# =======================================================================================================
#     Index alignment BAM files
# =======================================================================================================

rule align_reads_bam_index:
    input:
        bam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam",
        bam_sort_done = "results/align_reads_with_kmers/{phenos_filt}/bam_sort.done"
    output:
        bai = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam.bai",
        bam_index_done = touch("results/align_reads_with_kmers/{phenos_filt}/bam_index.done")
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
        samtools index -@ {threads} {input.bam} 2> {log}
        """

# =======================================================================================================
#     Convert alignment BAM files to BED format
# =======================================================================================================

rule align_reads_bam_to_bed:
    input:
        bam = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam",
        bai = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam.bai",
        bam_index_done = "results/align_reads_with_kmers/{phenos_filt}/bam_index.done"
    output:
        bed = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bed",
        bam_to_bed_done = touch("results/align_reads_with_kmers/{phenos_filt}/bam_to_bed.done")
    conda:
        "../envs/bedtools.yaml"
    threads:
        config["params"]["samtools"]["threads"]
    log:
        "logs/align_reads/{phenos_filt}/align_reads.bam_to_bed.log"
    message:
        "Converting alignment BAM files to BED format..."
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output.bed} 2> {log}
        """

# =======================================================================================================
#     Merge overlapping features in alignment BED files
# =======================================================================================================

rule merge_align_reads_bed:
    input:
        bed = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bed",
        bam_to_bed_done = "results/align_reads_with_kmers/{phenos_filt}/bam_to_bed.done"
    output:
        merged_bed = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.merged.bed",
        merge_bed_done = touch("results/align_reads_with_kmers/{phenos_filt}/merge_bed.done")
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/align_reads/{phenos_filt}/align_reads.merge_bed.log"
    message:
        "Merging alignment BED files..."
    shell:
        """
        bedtools merge -i {input.bed} -s -c 6 -o distinct > {output.merged_bed} 2> {log}
        """

# =======================================================================================================
#     tabix BED files
# =======================================================================================================

rule tabix_align_reads_bed:
    input:
        bed = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.merged.bed",
        merge_bed_done = "results/align_reads_with_kmers/{phenos_filt}/merge_bed.done"
    output:
        tbi = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.merged.bed.gz.tbi",
        tabix_bed_done = touch("results/align_reads_with_kmers/{phenos_filt}/tabix_bed.done")
    conda:
        "../envs/align_reads.yaml"
    log:
        "logs/align_reads/{phenos_filt}/align_reads.tabix_bed.log"
    message:
        "tabix BED files..."
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | bgzip > {input.bed}.gz 2> {log}
        tabix -p bed {input.bed}.gz 2> {log}
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