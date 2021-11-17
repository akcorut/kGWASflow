import csv
import os, argparse, glob, shutil
import pandas as pd
import subprocess 
import sys

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

        ##############################################################
        ############ Align reads with significant k-mers #############
        ##############################################################

    in_path = snakemake.params.in_prefix # Input path for reads with k-mers
    pheno = snakemake.params.pheno # Pnenotype name

    # Fetching paths for reads (R1 and R2) with k-mers for a given phenotype
    R1_list = glob.glob(in_path + '/' + pheno + '/' + "*_R1.fastq")

    R2_list = glob.glob(in_path + '/' + pheno + '/' + "*_R2.fastq")

    # Generating list of samples that have reads with k-mers for a given phenotype
    accession_list=[]
    for fastq in R1_list:
        fname = os.path.basename(fastq)
        accession_list.append(fname.split('_reads_with_kmers')[0])

    # Aligning reads for each sample
    for acc in accession_list:
        r1_path= in_path + '/' + pheno + '/' + acc + "_reads_with_kmers_R1.fastq"
        r2_path= in_path + '/' + pheno + '/' + acc + "_reads_with_kmers_R2.fastq"

        # Making output directory if doesn't exist
        if not os.path.exists(snakemake.params.out_prefix + '/' + pheno):
            os.makedirs(snakemake.params.out_prefix + '/' + pheno)

        print("################## Aligning reads with k-mers using bowtie2... ##################")
        print("")
        print("Sample: " + acc)
        print("")
        print("R1: " + r1_path)
        print("R2: " + r2_path)
        print("")
        print("")
        subprocess.run(" bowtie2 -p {threads} --very-sensitive-local -x {index} -1 {r1} -2 {r2} -S {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.sam 2> {log}/{accession}.reads_with_kmers_alignment.bowtie2.log".format(
            r1= r1_path, r2= r2_path, index= snakemake.params.index, 
            pheno=pheno, outdir=snakemake.params.out_prefix, 
            accession= acc, threads=snakemake.threads, log=snakemake.params.bowtie_log), shell=True)

        print("################## Converting SAM to BAM... ##################")
        print("")
        print("Sample: " + acc)
        print("")
        subprocess.run(" samtools view -S -b {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.sam > {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.bam".format(
            pheno=pheno, outdir=snakemake.params.out_prefix, accession= acc), shell=True)
        print("")

        print("################## Sorting alignment BAM... ##################")
        print("")
        print("Sample: " + acc)
        print("")
        subprocess.run(" samtools sort {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.bam -o {outdir}/{pheno}/{accession}_reads_with_kmers_alignment_sorted.bam".format(
            pheno=pheno, outdir=snakemake.params.out_prefix, accession= acc), shell=True)
        print("")    

        print("################## Indexing alignment BAM... ##################")
        print("")
        print("Sample: " + acc)
        print("")
        subprocess.run(" samtools index {outdir}/{pheno}/{accession}_reads_with_kmers_alignment_sorted.bam".format(
            pheno=pheno, outdir=snakemake.params.out_prefix, accession= acc), shell=True)
        print("")

        print("################## Removing alignment SAM... ##################")
        print("")
        print("Sample: " + acc)
        print("")   
        subprocess.run(" rm {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.sam".format(
            pheno=pheno, outdir=snakemake.params.out_prefix, accession= acc), shell=True)

        subprocess.run(" rm {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.bam".format(
            pheno=pheno, outdir=snakemake.params.out_prefix, accession= acc), shell=True)
        print("")
        print("Aligning reads with k-mers is done for the sample: " + acc)
        print("")

    print("")
    print("All alignments DONE.")