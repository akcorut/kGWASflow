# Configuration Settings

To configure this workflow, you need to edit `config/config.yaml` file based on your needs. The explanation for each parameter is provided in the `config.yaml` file. Below parameters has to be set in order to run this workflow:

- `samples`: the path of the sample sheet
- `phenotype`: the path of the phenotype sheet
- `ref:` `fasta:` the path of the refrence genome fasta file (Required for alignment step)

# Sample Sheet (samples.tsv)

Add samples to `config/samples.tsv`. For each sample, you need to specify the library name under `library_name` column. You can add one or more libraries for each sample. You also need to add either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ file path (can be any path in your csystem) for each sample-library combination.

**Note:** `SRA` option currently is not available therefore it can be left blank.


# Phenotype Sheet (phenos.tsv)

Add phenotype information to `config/phenos.tsv`. For each phenotype you need to specify the name of the phenotype (column `pheno_name`) and the path of the  correspondig phenotype file (column `pheno_path`). 

The format of the phenotye file has to match the description in [the kmersGWAS manual](https://github.com/voichek/kmersGWAS/blob/master/manual.pdf):

> The phenotype file should be with two columns separated by tabs (\t), with “accession_id” and “phenotype_value” in the header. The first column lists the name of the individuals, as in the kmers_table, and the second, phenotype_value, the phenotypic values in numerical form.
