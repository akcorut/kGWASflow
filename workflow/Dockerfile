FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="6cb40b386ca06bf0ba52e4db445cb08daba5c11395b28ba52c16e7ba6042fc6b"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/align_contigs.yaml
#   prefix: /conda-envs/de8c5ab797d5257f74c2e0ff78714b34
#   name: align_contigs
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - minimap2=2.22
#     - samtools=1.16.1
RUN mkdir -p /conda-envs/de8c5ab797d5257f74c2e0ff78714b34
COPY envs/align_contigs.yaml /conda-envs/de8c5ab797d5257f74c2e0ff78714b34/environment.yaml

# Conda environment:
#   source: envs/align_kmers.yaml
#   prefix: /conda-envs/862d702752dd97c0a8bceb0edb8b6570
#   name: align_kmers
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bowtie=1.3.1
#     - bowtie2=2.4.5
#     - samtools=1.16.1
RUN mkdir -p /conda-envs/862d702752dd97c0a8bceb0edb8b6570
COPY envs/align_kmers.yaml /conda-envs/862d702752dd97c0a8bceb0edb8b6570/environment.yaml

# Conda environment:
#   source: envs/align_reads.yaml
#   prefix: /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9
#   name: align_reads
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bowtie2=2.4.5
#     - samtools=1.16.1
#     - seqkit=0.16.1
RUN mkdir -p /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9
COPY envs/align_reads.yaml /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9/environment.yaml

# Conda environment:
#   source: envs/assemble_reads.yaml
#   prefix: /conda-envs/16900f85f28b56de8c694115f6622fdd
#   name: assemble_reads
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - spades=3.15
RUN mkdir -p /conda-envs/16900f85f28b56de8c694115f6622fdd
COPY envs/assemble_reads.yaml /conda-envs/16900f85f28b56de8c694115f6622fdd/environment.yaml

# Conda environment:
#   source: envs/build_kmers_gwas.yaml
#   prefix: /conda-envs/2aca5d40efbf8774681db1f06ade387b
#   name: build_kmers_gwas
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - wget=1.20.1
#     - unzip=6.0
RUN mkdir -p /conda-envs/2aca5d40efbf8774681db1f06ade387b
COPY envs/build_kmers_gwas.yaml /conda-envs/2aca5d40efbf8774681db1f06ade387b/environment.yaml

# Conda environment:
#   source: envs/kmc.yaml
#   prefix: /conda-envs/097d216c7fe7a0161751156cf35c077b
#   name: kmc
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - kmc=3.2.1
RUN mkdir -p /conda-envs/097d216c7fe7a0161751156cf35c077b
COPY envs/kmc.yaml /conda-envs/097d216c7fe7a0161751156cf35c077b/environment.yaml

# Conda environment:
#   source: envs/kmers_gwas.yaml
#   prefix: /conda-envs/ddd750164f6ced3aadf9f30d44f7d999
#   name: kmers_gwas
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - libgcc=7.2.0
RUN mkdir -p /conda-envs/ddd750164f6ced3aadf9f30d44f7d999
COPY envs/kmers_gwas.yaml /conda-envs/ddd750164f6ced3aadf9f30d44f7d999/environment.yaml

# Conda environment:
#   source: envs/kmers_gwas_py2.yaml
#   prefix: /conda-envs/b37c5dd23141def831b81ad5068da1d3
#   name: kmers_gwas_py2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - libiconv=1.17
#     - libgcc=7.2.0
#     - python=2.7.15
#     - r=3.5.1
#     - r-base=3.5.1
#     - r-mass=7.3
#     - r-matrixcalc=1.0
RUN mkdir -p /conda-envs/b37c5dd23141def831b81ad5068da1d3
COPY envs/kmers_gwas_py2.yaml /conda-envs/b37c5dd23141def831b81ad5068da1d3/environment.yaml

# Conda environment:
#   source: envs/kmers_stats.yaml
#   prefix: /conda-envs/ab577270450d84b1a5172d16e92a5d34
#   name: kmers_stats
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - matplotlib-base=3.2.2
#     - numpy=1.20.2
#     - python=3.7.0
#     - scipy=1.6.2
#     - seaborn=0.11.1
#     - seaborn-base=0.11.1
RUN mkdir -p /conda-envs/ab577270450d84b1a5172d16e92a5d34
COPY envs/kmers_stats.yaml /conda-envs/ab577270450d84b1a5172d16e92a5d34/environment.yaml

# Conda environment:
#   source: envs/plot_manhattan.yaml
#   prefix: /conda-envs/78847d20ad28fb2c74c621011802b919
#   name: plot_manhattan
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - matplotlib=3.6.2
#     - numpy=1.23.3
#     - scipy=1.9.3
#     - pandas=1.4.4
#     - seaborn=0.12.1
#     - natsort=7.1.1
#     - pip=22.2.2
#     - pip:
#       - qmplot
RUN mkdir -p /conda-envs/78847d20ad28fb2c74c621011802b919
COPY envs/plot_manhattan.yaml /conda-envs/78847d20ad28fb2c74c621011802b919/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/0.80.0/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/16114a9972039cb62b9656ba636b27b0
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bowtie2 ==2.4.4  # Keep consistent with version specified in bowtie2/align
#     - samtools ==1.10
RUN mkdir -p /conda-envs/16114a9972039cb62b9656ba636b27b0
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.80.0/bio/bowtie2/build/environment.yaml /conda-envs/16114a9972039cb62b9656ba636b27b0/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.12.2/bio/fastqc/environment.yaml
#   prefix: /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc ==0.11.9
RUN mkdir -p /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.12.2/bio/fastqc/environment.yaml /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.12.2/bio/multiqc/environment.yaml
#   prefix: /conda-envs/fd1008dfa88a4500724e3596f62f8bff
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.12
RUN mkdir -p /conda-envs/fd1008dfa88a4500724e3596f62f8bff
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.12.2/bio/multiqc/environment.yaml /conda-envs/fd1008dfa88a4500724e3596f62f8bff/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.12.2/bio/sra-tools/fasterq-dump/environment.yaml
#   prefix: /conda-envs/85633ff8bea713d372cb9152f291c3a8
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - sra-tools >2.9.1
#     - pigz >=2.6
#     - pbzip2 >=1.1
#     - snakemake-wrapper-utils =0.3
RUN mkdir -p /conda-envs/85633ff8bea713d372cb9152f291c3a8
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.12.2/bio/sra-tools/fasterq-dump/environment.yaml /conda-envs/85633ff8bea713d372cb9152f291c3a8/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/de8c5ab797d5257f74c2e0ff78714b34 --file /conda-envs/de8c5ab797d5257f74c2e0ff78714b34/environment.yaml && \
    mamba env create --prefix /conda-envs/862d702752dd97c0a8bceb0edb8b6570 --file /conda-envs/862d702752dd97c0a8bceb0edb8b6570/environment.yaml && \
    mamba env create --prefix /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9 --file /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9/environment.yaml && \
    mamba env create --prefix /conda-envs/16900f85f28b56de8c694115f6622fdd --file /conda-envs/16900f85f28b56de8c694115f6622fdd/environment.yaml && \
    mamba env create --prefix /conda-envs/2aca5d40efbf8774681db1f06ade387b --file /conda-envs/2aca5d40efbf8774681db1f06ade387b/environment.yaml && \
    mamba env create --prefix /conda-envs/097d216c7fe7a0161751156cf35c077b --file /conda-envs/097d216c7fe7a0161751156cf35c077b/environment.yaml && \
    mamba env create --prefix /conda-envs/ddd750164f6ced3aadf9f30d44f7d999 --file /conda-envs/ddd750164f6ced3aadf9f30d44f7d999/environment.yaml && \
    mamba env create --prefix /conda-envs/b37c5dd23141def831b81ad5068da1d3 --file /conda-envs/b37c5dd23141def831b81ad5068da1d3/environment.yaml && \
    mamba env create --prefix /conda-envs/ab577270450d84b1a5172d16e92a5d34 --file /conda-envs/ab577270450d84b1a5172d16e92a5d34/environment.yaml && \
    mamba env create --prefix /conda-envs/78847d20ad28fb2c74c621011802b919 --file /conda-envs/78847d20ad28fb2c74c621011802b919/environment.yaml && \
    mamba env create --prefix /conda-envs/16114a9972039cb62b9656ba636b27b0 --file /conda-envs/16114a9972039cb62b9656ba636b27b0/environment.yaml && \
    mamba env create --prefix /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8 --file /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8/environment.yaml && \
    mamba env create --prefix /conda-envs/fd1008dfa88a4500724e3596f62f8bff --file /conda-envs/fd1008dfa88a4500724e3596f62f8bff/environment.yaml && \
    mamba env create --prefix /conda-envs/85633ff8bea713d372cb9152f291c3a8 --file /conda-envs/85633ff8bea713d372cb9152f291c3a8/environment.yaml && \
    mamba clean --all -y
