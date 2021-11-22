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

    kmers_tab = snakemake.params.kmers_tab_prefix
    kmers_list = snakemake.params.kmers_list_prefix
    pheno = snakemake.params.pheno
    samp_tab = snakemake.input.samp_tab
    fetch_reads_dir = snakemake.params.fetch_reads_prefix
    

            ##############################################################
            ####### Identfiy accessions having the k-mers present ########
            ##############################################################

    ## Read "kmers_table.txt" for a given phenotype
    kmers_table = pd.read_csv(kmers_tab + "/" + pheno + "_kmers_table.txt", 
                            sep="\t")

    ## Filter out accessions/samples that doesn't have the kmers in them
    kmers_table_filt = kmers_table.loc[:, (kmers_table != 0).any(axis=0)]

    ## Get accesions with k-mers
    acc_list = list(kmers_table_filt.iloc[: , 1:].columns.values)

    print("Phenotype: " + pheno)
    print("")
    print("################## Identifying accessions that has the k-mers... ##################")
    print("")
    print(str(len(acc_list)) + " accessions identified: " + ', '.join(acc_list))

    print("")

    print("######################### Fetching reads with k-mers... #########################")

    ## Read samples sheet (samples.tsv)
    samples_tab = pd.read_csv(samp_tab, 
                                    sep="\t",
                                    dtype='object')

    ## Filter samples table to get only samples with k-mers
    samples_tab_filt = samples_tab[samples_tab['sample_name'].isin(acc_list)]
    samples_tab_filt = samples_tab_filt.reset_index(drop=True)

    ## Define the path kmers_list.fa file for given phenotype
    kmers_list_fa = kmers_list + "/" + pheno + "_kmers_list.fa"

    ## Check if output directory is exist. If not create one
    if not os.path.exists(args.out_dir + '/' + pheno):
        os.makedirs(args.out_dir + '/' + pheno)
    
    ## Filter the reads that have one of the k-mers in them
    for index, row in samples_tab_filt.iterrows():
        if row["sample_name"] == row["library_name"]:

            ## Source code: https://github.com/voichek/fetch_reads_with_kmers
            subprocess.run(" {fetch_reads_path}/fetch_reads {r1} {r2} {fa} 31 {out}/{pheno}/{acc}_reads_with_kmers".format(
                fetch_reads_path= fetch_reads_dir, out=args.out_dir, r1= row["fq1"], r2= row["fq2"],
                acc= row["sample_name"], lib = row["library_name"], fa= kmers_list_fa, pheno= pheno), shell=True)
        else:
            ## Source code: https://github.com/voichek/fetch_reads_with_kmers
            subprocess.run(" {fetch_reads_path}/fetch_reads {r1} {r2} {fa} 31 {out}/{pheno}/{acc}_{lib}_reads_with_kmers".format(
                fetch_reads_path= fetch_reads_dir, out=args.out_dir, r1= row["fq1"], r2= row["fq2"],
                acc= row["sample_name"], lib = row["library_name"], fa= kmers_list_fa, pheno= pheno), shell=True)