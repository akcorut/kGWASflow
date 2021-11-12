import csv
import sys
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description='')
parser.add_argument('-k', '--kmers_tab', help='Path for the filter_kmers results', type=str, required=True)
parser.add_argument('-l', '--kmers_list', help='Path for the fetch_kmers results', type=str, required=True)
parser.add_argument('-p', '--pheno', help='Phenotype name', type=str, required=True)
parser.add_argument('-s', '--samp_tab', help='Path for the samples sheet (samples.tsv)', type=str, required=True)
parser.add_argument('-o', '--out_dir', help='Output directory path', type=str, required=False)
parser.add_argument('-path', '--fetch_reads_path', help='Path for the fetch_reads_with_kmers library', type=str, required=True)
args = parser.parse_args()

        ##############################################################
        ####### Identfiy accessions having the k-mers present ########
        ##############################################################

## Read "kmers_table.txt" for a given phenotype
kmers_table = pd.read_csv(args.kmers_tab + "/" + args.pheno + "_kmers_table.txt", 
                        sep="\t")

## Filter out accessions that doesn't have the kmers in them
kmers_table_filt = kmers_table.loc[:, (kmers_table != 0).any(axis=0)]
# print(kmers_table_filt)

acc_list = list(kmers_table_filt.iloc[: , 1:].columns.values)

print("Phenotype: " + args.pheno)
print("")
print("################## Identifying accessions that has the k-mers... ##################")
print("")
print(str(len(acc_list)) + " accessions identified: " + ', '.join(acc_list))

print("")

print("######################### Fetching reads with k-mers... #########################")

samples_tab = pd.read_csv(args.samp_tab, 
                                sep="\t",
                                dtype='object')

samples_tab_filt = samples_tab[samples_tab['sample_name'].isin(acc_list)]
samples_tab_filt = samples_tab_filt.reset_index(drop=True)

## Define the path kmers_list.fa file for given phenotype
kmers_list_fa = args.kmers_list + "/" + args.pheno + "_kmers_list.fa"

## Check if output directory is exist. If not create one
if not os.path.exists(args.out_dir + '/' + args.pheno):
    os.makedirs(args.out_dir + '/' + args.pheno)

# for index, row in samples_tab_filt.iterrows():
#         print(" fetch_reads {r1} {r2} 31 {acc}_{lib}_reads_with_kmers".format( 
#                 r1=row["fq1"], r2=row["fq2"], acc=row["accession"], lib=row["library"]))
        
## Filter the reads that have one of the k-mers in them
for index, row in samples_tab_filt.iterrows():
    if row["sample_name"] == row["library_name"]:
        subprocess.run(" {fetch_reads_path}/fetch_reads {r1} {r2} {fa} 31 {out}/{pheno}/{acc}_reads_with_kmers".format(
            fetch_reads_path= args.fetch_reads_path, out=args.out_dir, r1= row["fq1"], r2= row["fq2"],
            acc= row["sample_name"], lib = row["library_name"], fa= kmers_list_fa, pheno= args.pheno), shell=True)
    else:
        subprocess.run(" {fetch_reads_path}/fetch_reads {r1} {r2} {fa} 31 {out}/{pheno}/{acc}_{lib}_reads_with_kmers".format(
            fetch_reads_path= args.fetch_reads_path, out=args.out_dir, r1= row["fq1"], r2= row["fq2"],
            acc= row["sample_name"], lib = row["library_name"], fa= kmers_list_fa, pheno= args.pheno), shell=True)