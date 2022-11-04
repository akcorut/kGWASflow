# Configuration Settings

To configure this workflow, you need to edit `config/config.yaml` file based on your needs. The explanation for each parameter is provided in the `config.yaml` file. Below parameters has to be set in order to run this workflow:

- `samples`: the path of the sample sheet
- `phenotype`: the path of the phenotype sheet
- `ref:` `fasta:` the path of the refrence genome fasta file (Required for Post-GWAS mapping steps)

# Sample Sheet (samples.tsv)

Add samples information to `config/samples.tsv`. For each sample, you need to specify the sample name under the `sample_name` column and the corresponding library name under the `library_name` column.

* Each sample has a `library_name`, which can be e.g. a library name, sequencing run or lane.
* Each sample has a `sample_name`, which associates it with the biological sample it comes from.
* For each sample, define the paths for your (columns `fq1`, `fq2`) FASTQ files (these can point to any path in your system). 
* As an alternative, you can provide an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. If ypu provde SRA accession instead of FASTQ file paths, the workflow will automatically download the corresponding paired end reads from SRA. If both local file paths and SRA accession are provided, the local files will be preferred.

## Example samples.tsv with local files:

| sample_name | library_name | fq1 | fq2 | sra |
|---|---|---|---|---|
| 33-16 | HB1JEADXX_GATCAG | /path/to/your/local/fastq_1/33-16_HB1JEADXX_GATCAG_1.fastq.gz | /path/to/your/local/fastq_2/33-16_HB1JEADXX_GATCAG_2.fastq.gz |  |
| 33-16 | HB1LCADXX_GATCAG | /path/to/your/local/fastq_1/33-16_HB1LCADXX_GATCAG_1.fastq.gz | /path/to/your/local/fastq_2/33-16_HB1LCADXX_GATCAG_2.fastq.gz |  |
| 33-16 | HLC27CCXX_L7 | /path/to/your/local/fastq_1/33-16/33-16_HLC27CCXX_L7_1.fastq.gz | /path/to/your/local/fastq_2/33-16_HLC27CCXX_L7_2.fastq.gz |  |
| 38-11 | HLCY2CCXX_L3 | /path/to/your/local/fastq_1/38-11/38-11_HLCY2CCXX_L3_1.fastq.gz | /path/to/your/local/fastq_2/38-11_HLCY2CCXX_L3_2.fastq.gz |  |
| 4226 | H174VBGXX | /path/to/your/local/fastq_1/4226/4226_H174VBGXX_1.fastq.gz | /path/to/your/local/fastq_2/4226_H174VBGXX_2.fastq.gz |  |
| 4226 | HLC27CCXX_L3 | /path/to/your/local/fastq_1/4226/4226_HLC27CCXX_L3_1.fastq.gz | /path/to/your/local/fastq_2/4226/4226_HLC27CCXX_L3_2.fastq.gz |  |
| 4722 | H3CWKBGXX_GCCAAT | /path/to/your/local/fastq_1/4722_H3CWKBGXX_GCCAAT_1.fastq.gz | /path/to/your/local/fastq_2/4722_H3CWKBGXX_GCCAAT_2.fastq.gz |  |
| 4722 | HLCY2CCXX_L3 | /path/to/your/local/fastq_1/4722_HLCY2CCXX_L3_1.fastq.gz | /path/to/your/local/fastq_2/4722_HLCY2CCXX_L3_2.fastq.gz |  |
| A188 | HLC27CCXX_L3 | /path/to/your/local/fastq_1/A188_HLC27CCXX_L3_1.fastq.gz | /path/to/your/local/fastq_2/A188_HLC27CCXX_L3_2.fastq.gz |  |
| A188 | H174VBGXX | /path/to/your/local/fastq_1/A188_H174VBGXX_1.fastq.gz | /path/to/your/local/fastq_2/A188_H174VBGXX_2.fastq.gz |  |

## Example samples.tsv with SRA accessions:

| sample_name | library_name | fq1 | fq2 | sra |
|---|---|---|---|---|
| C00022105 | C00022105 |  |  | SRR3050845 |
| C00022106 | C00022106 |  |  | SRR3050846 |
| C00022107 | C00022107 |  |  | SRR3050847 |
| C00022108 | C00022108 |  |  | SRR3050848 |
| C00022109 | C00022109 |  |  | SRR3050849 |
| C00022110 | C00022110 |  |  | SRR3050850 |
| C00022119 | C00022119 |  |  | SRR3050851 |
| C00022112 | C00022112 |  |  | SRR3050852 |
| C00022113 | C00022113 |  |  | SRR3050853 |
| C00022114 | C00022114 |  |  | SRR3050854 |

# Phenotype Sheet (phenos.tsv)

Add phenotype information to `config/phenos.tsv`. For each phenotype you need to specify the name of the phenotype (column `pheno_name`) and the path of the  correspondig phenotype file (column `pheno_path`). 

## Example phenos.tsv:

Below is an example **phenos.tsv** file. First columns corresponds to phenotype name and the second column corresponds to its file path.

| pheno_name | pheno_path |
|---|---|
| Flowering_Time | path/to/pheno/Flowering_Time.pheno |
| Seed_Oil | path/to/pheno/Seed_Oil.pheno |
| UpperLeaf_Angle | path/to/pheno/UpperLeaf_Angle.pheno |

**Note:** The name in in the pheno_name column must match to name of the .pheno file. For example, if pheno_name name is "Flowering_Time", .phone file must be named as Flowering_Time.pheno. See above table for examples.

## Example Phenotype File (pheno_name.pheno)

The format of the phenotye file has to match the description in [the kmersGWAS manual](https://github.com/voichek/kmersGWAS/blob/master/manual.pdf):

> The phenotype file should be with two columns separated by tabs (\t), with “accession_id” and “phenotype_value” in the header. The first column lists the name of the individuals, as in the kmers_table, and the second, phenotype_value, the phenotypic values in numerical form.

Below is an example of a **phenotype file** (pheno_name.pheno):

| accession_id | phenotype_value |
|---|---|
| 33-16 | 81.704 |
| 38-11 | 83.33 |
| 4226 | 73.066 |
| 4722 | 83.016 |
| A188 | 67.996 |
| A214N | 87.855 |
| A239 | 73.357 |

