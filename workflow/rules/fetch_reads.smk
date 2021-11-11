# =======================================================================================================
#     Download and compile fetch_reads_with_kmers 
#     (Source: https://github.com/voichek/fetch_reads_with_kmers)
# =======================================================================================================


rule download_fetch_reads_with_kmers:
    output:
        cpp= "scripts/external/fetch_reads_with_kmers/fetch_reads.cpp",
    shell:
        """
        git clone --recurse-submodules https://github.com/voichek/fetch_reads_with_kmers.git
        mv fetch_reads_with_kmers /scripts/external
        """

rule make_fetch_reads_with_kmers:
    input:
        rules.download_fetch_reads_with_kmers.output.cpp
    output:
        "scripts/external/fetch_reads_with_kmers/fetch_reads"
    params:
        prefix= basedir + "/scripts/fetch_reads_with_kmers"
    shell:
        """
        cd {params.prefix}
        make
        """