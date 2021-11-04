# =================================================================================================
#     Download kmersGWAS
# =================================================================================================

rule download_kmersGWAS:
    output:
        temp(KMERSGWAS_ZIP_PATH)
    conda:
        "../envs/build_kmers_gwas.yaml"
    message: 
        "Downloading kmersGWAS source code..."
    shell:
        "wget https://github.com/voichek/kmersGWAS/releases/download/{KMERSGWAS_VERSION}/{KMERSGWAS_ZIP_PREFIX}.zip -O {output}"

# =================================================================================================
#     Extract kmersGWAS
# =================================================================================================

rule extract_kmersGWAS:
    input:
        rules.download_kmersGWAS.output
    output:
        kmersGWAS_py = KMERSGWAS_PY_PATH,
        kmersGWAS_bin = directory(KMERSGWAS_BIN_PATH)
    conda:
        "../envs/build_kmers_gwas.yaml"
    message: 
        "Unzipping kmersGWAS source code..."
    shell:
        "unzip {input} -d {KMERSGWAS_DIR}"

# =================================================================================================