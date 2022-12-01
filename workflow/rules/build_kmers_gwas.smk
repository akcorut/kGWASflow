# =================================================================================================
#     Download kmersGWAS
# =================================================================================================

rule download_kmersGWAS:
    output:
        temp(KMERSGWAS_ZIP_PATH)
    params:
        version = KMERSGWAS_VERSION,
        prefix= KMERSGWAS_ZIP_PREFIX
    conda:
        "../envs/build_kmers_gwas.yaml"
    log:
        "logs/build_kmers_gwas/download_kmersGWAS.wget.log"
    message: 
        "Downloading kmersGWAS source code..."
    shell:
        "wget https://github.com/voichek/kmersGWAS/releases/download/{params.version}/{params.prefix}.zip -O {output} 2> {log}"

# =================================================================================================
#     Extract kmersGWAS
# =================================================================================================

rule extract_kmersGWAS:
    input:
        rules.download_kmersGWAS.output
    output:
        kmersGWAS_py = KMERSGWAS_PY_PATH,
        kmersGWAS_bin = directory(KMERSGWAS_BIN_PATH)
    params:
        dir= lambda w, input: os.path.dirname(input[0])
    conda:
        "../envs/build_kmers_gwas.yaml"
    log:
        "logs/build_kmers_gwas/extract_kmersGWAS.unzip.log"
    message: 
        "Unzipping kmersGWAS source code..."
    shell:
        "unzip {input} -d {params.dir} 2> {log}"

# =================================================================================================