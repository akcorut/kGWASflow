# =================================================================================================
#     Generate IGV HTML report for contigs alignment
# =================================================================================================

rule align_contigs_igv_report:
    input:
        fasta="resources/ref/genome/genome.fasta",
        fai="resources/ref/genome/genome.fasta.fai",
        vcf="results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bed",
        tracks= [
            get_genome_annotation(),
            "results/align_contigs/{phenos_filt}/alignment/{phenos_filt}_contigs_aligned.filter.sorted.bam",
        ],
    output:
        igv_report= report(
            "results/igv_reports/align_contigs/{phenos_filt}/{phenos_filt}_contigs_aligned.igv-report.html",
            caption= "../report/generate_igv_report_contigs.rst",
            category= "IGV Reports - Contig Mapping",
            subcategory="{phenos_filt}",
        )
    params:
        extra= config["params"]["igv_report"]["extra"] # optional params
    log:
        "logs/igv_report/align_contigs/{phenos_filt}.igv-report.log"
    wrapper:
        "v1.25.0/bio/igv-reports"

# =========================================================================================================
#     Aggregate igv_report outputs
# =========================================================================================================

rule aggregate_align_contigs_igv_reports:
    input:
        aggregate_input_align_contigs_igv_report
    output:
        "results/igv_reports/align_contigs.igv_report.done"
    log:
        "logs/igv_report/aggregate_align_contigs_igv_reports.log"
    message:
        "Checking if generating igv_report is done..."
    shell:
        """
        touch {output} 2> {log}
        """
    
# =================================================================================================
#     Generate IGV HTML report for reads alignments
# =================================================================================================

rule align_reads_igv_report:
    input:
        fasta="resources/ref/genome/genome.fasta",
        fai="resources/ref/genome/genome.fasta.fai",
        vcf="results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.merged.bed",
        tbi = "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.merged.bed.gz.tbi",
        tracks= [
            get_genome_annotation(),
            "results/align_reads_with_kmers/{phenos_filt}/{phenos_filt}.align_reads_with_kmers.filter.sorted.bam",
        ],
        done = "results/align_reads_with_kmers/{phenos_filt}/bam_to_bed.done"
    output:
        igv_report= report(
            "results/igv_reports/align_reads/{phenos_filt}/{phenos_filt}_reads_aligned.igv-report.html",
            caption= "../report/generate_igv_report_reads.rst",
            category= "IGV Reports - Read Mapping",
            subcategory="{phenos_filt}",
        )
    params:
        extra= config["params"]["igv_report"]["extra"] # optional params
    threads:
        config["params"]["igv_report"]["threads"]
    log:
        "logs/igv_report/align_reads/{phenos_filt}.igv-report.log"
    wrapper:
        "v1.25.0/bio/igv-reports"

# =========================================================================================================
#     Aggregate igv_report outputs
# =========================================================================================================

rule aggregate_align_reads_igv_reports:
    input:
        aggregate_input_align_reads_igv_report
    output:
        "results/igv_reports/align_reads.igv_report.done"
    log:
        "logs/igv_report/aggregate_align_reads_igv_reports.log"
    message:
        "Checking if generating igv_report is done..."
    shell:
        """
        touch {output} 2> {log}
        """

# =========================================================================================================