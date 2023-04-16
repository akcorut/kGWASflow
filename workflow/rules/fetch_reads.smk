# =======================================================================================================
#     Download and compile fetch_reads_with_kmers 
#     # Source: https://github.com/voichek/fetch_reads_with_kmers
# =======================================================================================================

rule download_fetch_reads_with_kmers:
    output:
        cpp= "scripts/external/fetch_reads_with_kmers/fetch_reads.cpp",
    log:
        log1 = "logs/fetch_reads/clone_fetch_reads_with_kmers.log",
        log2 = "logs/fetch_reads/move_fetch_reads_with_kmers.log"
    shell:
        """
        git clone --recurse-submodules https://github.com/voichek/fetch_reads_with_kmers.git 2> {log.log1}
        mv fetch_reads_with_kmers scripts/external 2> {log.log2}
        """

rule make_fetch_reads_with_kmers:
    input:
        rules.download_fetch_reads_with_kmers.output.cpp
    output:
        "scripts/external/fetch_reads_with_kmers/fetch_reads"
    params:
        prefix= lambda w, output: os.path.dirname(output[0])
    log:
        log1 = "logs/fetch_reads/make_fetch_reads_with_kmers.log",
    shell:
        """
        cd {params.prefix}
        make
        """

# =======================================================================================================
#     Fetch source reads of significant k-mers
# =======================================================================================================

if not config["settings"]["trimming"]["activate"]:
    rule fetch_source_reads:
        input:
            kmers_tab = "results/filter_kmers/{phenos_filt}_kmers_table.txt",
            kmers_list = "results/fetch_kmers/{phenos_filt}_kmers_list.fa",
            fetch_reads = "scripts/external/fetch_reads_with_kmers/fetch_reads",
            filter_kmers = "results/filter_kmers/filter_kmers.done",
            r1 = expand(rules.create_symlink.output.r1, zip,
                            sample=sample_names,
                            library=library_names),
            r2 = expand(rules.create_symlink.output.r2, zip,
                            sample=sample_names,
                            library=library_names),
        output:
            out_dir = directory("results/fetch_reads_with_kmers/{phenos_filt}/individual_reads"),
            done = touch("results/fetch_reads_with_kmers/{phenos_filt}/individual_reads/{phenos_filt}.fetching_source_reads.done")
        params:
            kmers_list_prefix = lambda w, input: os.path.dirname(input.kmers_list),
            samples= sample_names,
            library= library_names,
            pheno = "{phenos_filt}",
            kmer_len = config["params"]["kmc"]["kmer_len"]
        log:
            "logs/fetch_reads/{phenos_filt}/fetch_source_reads_of_kmers.log"
        threads: 
            config["params"]["fetch_reads"]["threads"]
        message:
            "Fetching reads that contain significant k-mers find in {input.kmers_list}..."    
        script:
            "../scripts/fetch_source_reads.py"

if config["settings"]["trimming"]["activate"]:
    rule fetch_source_reads:
        input:
            kmers_tab = "results/filter_kmers/{phenos_filt}_kmers_table.txt",
            kmers_list = "results/fetch_kmers/{phenos_filt}_kmers_list.fa",
            fetch_reads = "scripts/external/fetch_reads_with_kmers/fetch_reads",
            filter_kmers = "results/filter_kmers/filter_kmers.done",
            r1 = expand(rules.cutadapt_pe.output.fastq1, zip,
                            sample=sample_names,
                            library=library_names),
            r2 = expand(rules.cutadapt_pe.output.fastq2, zip,
                            sample=sample_names,
                            library=library_names),
        output:
            out_dir = directory("results/fetch_reads_with_kmers/{phenos_filt}/individual_reads"),
            done = touch("results/fetch_reads_with_kmers/{phenos_filt}/individual_reads/{phenos_filt}.fetching_source_reads.done")
        params:
            kmers_list_prefix = lambda w, input: os.path.dirname(input.kmers_list),
            samples= sample_names,
            library= library_names,
            pheno = "{phenos_filt}",
            kmer_len = config["params"]["kmc"]["kmer_len"]
        log:
            "logs/fetch_reads/{phenos_filt}/fetch_source_reads_of_kmers.log"
        threads: 
            config["params"]["fetch_reads"]["threads"]
        message:
            "Fetching reads that contain significant k-mers find in {input.kmers_list}..."    
        script:
            "../scripts/fetch_source_reads.py"


# =========================================================================================================
#     Aggregate fetch_source_reads outputs
# =========================================================================================================

rule aggregate_fetch_reads:
    input:
        aggregate_input_fetch_reads
    output:
        "results/fetch_reads_with_kmers/fetch_source_reads.done"
    log:
        "logs/fetch_reads/aggregate_fetch_reads.log"
    message:
        "Checking if fetching source reads of k-mers is done..."
    shell:
        """
        touch {output} 2> {log}
        """

# =========================================================================================================