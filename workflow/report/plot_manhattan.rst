Manhattan plot showing significantly associated k-mers (only k-mers pass p-value threshold) and their genomic locations. k-mers were mapped to a reference genome using {% if snakemake.config["settings"]["align_kmers"]["use_bowtie"] %} bowtie_ {% endif %} {% if snakemake.config["settings"]["align_kmers"]["use_bowtie2"] %} bowtie2_ {% endif %}. Manhattan plots were generated using the python library QMplot_ .

.. _bowtie: https://bowtie-bio.sourceforge.net/index.shtml
.. _bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _QMplot: https://github.com/ShujiaHuang/qmplot