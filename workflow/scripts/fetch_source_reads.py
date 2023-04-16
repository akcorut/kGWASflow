import csv
import sys
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
                #################################################
                ####### Load Snakemake Inputs and Params ########
                #################################################
    
    fq1_list= list(snakemake.input.r1) # list of fastq r1 files
    fq2_list= list(snakemake.input.r2) # list of fastq r2 files
    samp_list = list(snakemake.params.samples) # list of sample names
    lib_list = list(snakemake.params.library) # list of library names
    
    kmers_tab = snakemake.input.kmers_tab # kmers table

    kmers_list_fa = snakemake.input.kmers_list # kmers list fasta

    pheno = snakemake.params.pheno # phenotype name

    out_dir = snakemake.output.out_dir # output directory

    kmer_len = snakemake.params.kmer_len # kmer length

    fetch_reads = snakemake.input.fetch_reads # fetch_reads script

        
                    #########################################
                    ####### Generate a Samples Table ########
                    #########################################
    
    ## Make a samples info table w/ sample names, library names and corresponding fastq file paths
    reads_path_tab = {'sample_name': samp_list, 'library_name': lib_list, 'fq1': fq1_list, 'fq2': fq2_list}
    reads_path_tab = pd.DataFrame(data=reads_path_tab, dtype=object)
    print(reads_path_tab.head())

    # reads_path_tab.to_csv(snakemake.params.out_prefix + "test.txt", index=False, sep="\t")

            ###############################################################
            ####### Identify the samples having the k-mers present ########
            ###############################################################

    ## Read "kmers_table.txt" for a given phenotype
    kmers_table = pd.read_csv(kmers_tab, 
                            sep="\t")

    ## Filter out accessions/samples that doesn't have the kmers in them
    kmers_table_filt = kmers_table.loc[:, (kmers_table != 0).any(axis=0)]

    ## Get accesions with k-mers
    acc_list = list(kmers_table_filt.iloc[: , 1:].columns.values)

    print("Working on the phenotype: " + pheno)
    print("")
    print("################## Identifying accessions that has the k-mers... ##################")
    print("")
    print(str(len(acc_list)) + " accessions identified: " + ', \n'.join(acc_list))
    print("")
    print("######################### Fetching reads with k-mers... #########################")
    print("")
    
    ## Filter samples table to get only samples with k-mers
    reads_path_tab_filt = reads_path_tab[reads_path_tab['sample_name'].isin(acc_list)]
    reads_path_tab_filt = reads_path_tab_filt.reset_index(drop=True)

    ## Check if output directory is exist. If not create one
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    ## Filter the reads that have one of the k-mers in them
    for index, row in reads_path_tab_filt.iterrows():
        if row["sample_name"] == row["library_name"]:

            ## Source code: https://github.com/voichek/fetch_reads_with_kmers
            subprocess.run(" {fetch_reads} {r1} {r2} {fa} {kmer_len} {out}/{acc}_reads_with_kmers".format(
                fetch_reads= fetch_reads, out=out_dir, 
                r1= row["fq1"], r2= row["fq2"],
                acc= row["sample_name"], lib = row["library_name"], 
                fa= kmers_list_fa, pheno= pheno, kmer_len= kmer_len), shell=True)
            
            # make sure that the output file is generated and if not exit the script with error
            if not os.path.exists("{out}/{acc}_reads_with_kmers_R1.fastq".format(out=out_dir, acc=row["sample_name"])):
                print("Error: {acc}_reads_with_kmers_R1.fastq is not generated!".format(acc=row["sample_name"]))
                sys.exit(1)
            if not os.path.exists("{out}/{acc}_reads_with_kmers_R2.fastq".format(out=out_dir, acc=row["sample_name"])):
                print("Error: {acc}_reads_with_kmers_R2.fastq is not generated!".format(acc=row["sample_name"]))
                sys.exit(1)
                
        else:
            ## Source code: https://github.com/voichek/fetch_reads_with_kmers
            subprocess.run(" {fetch_reads} {r1} {r2} {fa} {kmer_len} {out}/{acc}_{lib}_reads_with_kmers".format(
                fetch_reads= fetch_reads, out=out_dir, 
                r1= row["fq1"], r2= row["fq2"],
                acc= row["sample_name"], lib = row["library_name"], 
                fa= kmers_list_fa, pheno= pheno, kmer_len= kmer_len), shell=True)
            
             # make sure that the output file is generated and if not exit the script with error
            if not os.path.exists("{out}/{acc}_{lib}_reads_with_kmers_R1.fastq".format(out=out_dir, acc=row["sample_name"], lib=row["library_name"])): 
                print("Error: {acc}_{lib}_reads_with_kmers_R1.fastq is not generated!".format(acc=row["sample_name"], lib=row["library_name"]))
                sys.exit(1)
            if not os.path.exists("{out}/{acc}_{lib}_reads_with_kmers_R2.fastq".format(out=out_dir, acc=row["sample_name"], lib=row["library_name"])):
                print("Error: {acc}_{lib}_reads_with_kmers_R2.fastq is not generated!".format(acc=row["sample_name"], lib=row["library_name"]))
                sys.exit(1)
            
            
    