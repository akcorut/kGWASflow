FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="6efeb294fe994b2ceab92aba0d070b0787dd5383d1691bc6bb40d4f1d2974bef"

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
#   prefix: /conda-envs/92483642db5b370e34a4f85d1d6ae1ce
#   name: assemble_reads
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - spades=3.15.5
RUN mkdir -p /conda-envs/92483642db5b370e34a4f85d1d6ae1ce
COPY envs/assemble_reads.yaml /conda-envs/92483642db5b370e34a4f85d1d6ae1ce/environment.yaml

# Conda environment:
#   source: envs/bedtools.yaml
#   prefix: /conda-envs/8907b830c092d734a11f64bde6f31eb5
#   name: bedtools
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bedtools=2.30.0
RUN mkdir -p /conda-envs/8907b830c092d734a11f64bde6f31eb5
COPY envs/bedtools.yaml /conda-envs/8907b830c092d734a11f64bde6f31eb5/environment.yaml

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
#   source: envs/kmer_stats.yaml
#   prefix: /conda-envs/a624324b5522c2cac3507a99cd8ae550
#   name: kmers_stats
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
RUN mkdir -p /conda-envs/a624324b5522c2cac3507a99cd8ae550
COPY envs/kmer_stats.yaml /conda-envs/a624324b5522c2cac3507a99cd8ae550/environment.yaml

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
#   source: envs/plot_manhattan.yaml
#   prefix: /conda-envs/58aaa9fc267769136b0866fbe323ee3c
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
#     - pysam=0.20.0
#     - pip:
#       - qmplot
RUN mkdir -p /conda-envs/58aaa9fc267769136b0866fbe323ee3c
COPY envs/plot_manhattan.yaml /conda-envs/58aaa9fc267769136b0866fbe323ee3c/environment.yaml

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
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.23.5/bio/sra-tools/fasterq-dump/environment.yaml
#   prefix: /conda-envs/09551784d42806e996c779649a93ea18
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - sra-tools =3.0.3
#     - pigz =2.6
#     - pbzip2 =1.1.13
#     - snakemake-wrapper-utils =0.5.0
RUN mkdir -p /conda-envs/09551784d42806e996c779649a93ea18
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.23.5/bio/sra-tools/fasterq-dump/environment.yaml /conda-envs/09551784d42806e996c779649a93ea18/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.25.0/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/c0b26e5a843c8355257bd60e19bbbaff
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.1
RUN mkdir -p /conda-envs/c0b26e5a843c8355257bd60e19bbbaff
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.25.0/bio/bowtie2/build/environment.yaml /conda-envs/c0b26e5a843c8355257bd60e19bbbaff/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.25.0/bio/igv-reports/environment.yaml
#   prefix: /conda-envs/a265d01756075df3edd511d305122ae5
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - igv-reports =1.7.0
RUN mkdir -p /conda-envs/a265d01756075df3edd511d305122ae5
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.25.0/bio/igv-reports/environment.yaml /conda-envs/a265d01756075df3edd511d305122ae5/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/de8c5ab797d5257f74c2e0ff78714b34 --file /conda-envs/de8c5ab797d5257f74c2e0ff78714b34/environment.yaml && \
    mamba env create --prefix /conda-envs/862d702752dd97c0a8bceb0edb8b6570 --file /conda-envs/862d702752dd97c0a8bceb0edb8b6570/environment.yaml && \
    mamba env create --prefix /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9 --file /conda-envs/b1ead5c69bc2d9b74b6c3aa85466f7c9/environment.yaml && \
    mamba env create --prefix /conda-envs/92483642db5b370e34a4f85d1d6ae1ce --file /conda-envs/92483642db5b370e34a4f85d1d6ae1ce/environment.yaml && \
    mamba env create --prefix /conda-envs/8907b830c092d734a11f64bde6f31eb5 --file /conda-envs/8907b830c092d734a11f64bde6f31eb5/environment.yaml && \
    mamba env create --prefix /conda-envs/2aca5d40efbf8774681db1f06ade387b --file /conda-envs/2aca5d40efbf8774681db1f06ade387b/environment.yaml && \
    mamba env create --prefix /conda-envs/097d216c7fe7a0161751156cf35c077b --file /conda-envs/097d216c7fe7a0161751156cf35c077b/environment.yaml && \
    mamba env create --prefix /conda-envs/a624324b5522c2cac3507a99cd8ae550 --file /conda-envs/a624324b5522c2cac3507a99cd8ae550/environment.yaml && \
    mamba env create --prefix /conda-envs/ddd750164f6ced3aadf9f30d44f7d999 --file /conda-envs/ddd750164f6ced3aadf9f30d44f7d999/environment.yaml && \
    mamba env create --prefix /conda-envs/b37c5dd23141def831b81ad5068da1d3 --file /conda-envs/b37c5dd23141def831b81ad5068da1d3/environment.yaml && \
    mamba env create --prefix /conda-envs/58aaa9fc267769136b0866fbe323ee3c --file /conda-envs/58aaa9fc267769136b0866fbe323ee3c/environment.yaml && \
    mamba env create --prefix /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8 --file /conda-envs/573927d1a2f1de4bfdd03a5385f50ed8/environment.yaml && \
    mamba env create --prefix /conda-envs/fd1008dfa88a4500724e3596f62f8bff --file /conda-envs/fd1008dfa88a4500724e3596f62f8bff/environment.yaml && \
    mamba env create --prefix /conda-envs/09551784d42806e996c779649a93ea18 --file /conda-envs/09551784d42806e996c779649a93ea18/environment.yaml && \
    mamba env create --prefix /conda-envs/c0b26e5a843c8355257bd60e19bbbaff --file /conda-envs/c0b26e5a843c8355257bd60e19bbbaff/environment.yaml && \
    mamba env create --prefix /conda-envs/a265d01756075df3edd511d305122ae5 --file /conda-envs/a265d01756075df3edd511d305122ae5/environment.yaml && \
    mamba clean --all -y
