.. image:: https://user-images.githubusercontent.com/42179487/194161153-cc832e57-dd03-481b-8eed-34cb13ba3097.png
    :width: 400
    :align: center

**kGWASflow:  A Snakemake Pipeline for k-mers-based GWAS** 

This workflow performs *k*-mers-based GWAS analysis on paired-end sequencing data. Quality control analysis was done using `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ and `MultiQC <https://multiqc.info>`_. 

{% if snakemake.config["settings"]["trimming"]["activate"] %}
Read trimming were performed using `Cutadapt <http://cutadapt.readthedocs.io>`_.
{% endif %}

*k*-mers were counted using `KMC <https://github.com/refresh-bio/KMC>`_. *k*-mers length was set to k= {{ snakemake.config["params"]["kmc"]["kmer_len"] }}. *k*-mers were counted twice with KMC, once with canonization and once without canonization as described in `Voichek et al. (2020) <https://www.nature.com/articles/s41588-020-0612-7>`_. After *k*-mers counting step, the *k*-mers from each sample are combined into a single binary file and filtered. The *k*-mers table with the presence/absence information of each *k*-mer was generated and kmersGWAS was performed using `kmersGWAS library <https://github.com/voichek/kmersGWAS>`_. *k*-mers-based GWAS was done by following `the method of Voichek et al. (2020) <https://github.com/voichek/kmersGWAS/blob/master/manual.pdf>`_.

---------------------------------------------

If you use **kGWASflow**, please **cite**:

* Kivanc Corut. akcorut/kGWASflow: v1.2.0. (2023). https://doi.org/10.5281/zenodo.7290926

* Voichek, Y., Weigel, D. Identifying genetic variants underlying phenotypic variation in plants without complete genomes. Nat Genet 52, 534–540 (2020). https://doi.org/10.1038/s41588-020-0612-7