import csv
import math
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import sys

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--in_tab', help='Input kmersGWAS results summary table', 
                    type=str, required=True)
parser.add_argument('-p', '--phenos', help='Phenotype', 
                    type=str, required=True)
parser.add_argument('-o', '--out_path', help='Output directory path', 
                    type=str, required=True)
parser.add_argument('-t', '--threshold', help='kmersGWAS FWER threshold for significant kmers (5 percent or 10 percent)', 
                    type=int, required=True, default=5)
args = parser.parse_args()

## Read phenotypes with significank k-mers table
phenos_list = pd.read_csv(args.phenos, 
                        sep="\t")

for pheno in phenos_list["pheno"]:
    ## Read the input kmersGWAS results summary table
    results_table = pd.read_csv(args.in_tab, 
                        sep="\t")
    
    results_table["kmers_pass_threshold_{thr}per".format(thr=str(args.threshold))] = results_table["kmers_pass_threshold_{thr}per".format(thr=str(args.threshold))].apply(eval)
    
    # Check if the current phnotype has any significant kmers associated with it
    if results_table.loc[results_table['pheno'] == pheno, 'num_kmers_pass_threshold_{thr}per'.format(thr=str(args.threshold))].iloc[0] == 0:
        pass
        print("No k-mers passing {thr}per threshold was found for the phenotype: ".format(thr=str(args.threshold)) + pheno)
    else:
        # If significant kmers are found, convert them into a list and write out the output as a table
        kmers_pass_threshold = results_table.loc[results_table['pheno'] == pheno, 'kmers_pass_threshold_{thr}per'.format(thr=str(args.threshold))]
    
        kmers_pass_threshold = kmers_pass_threshold.explode()

        print("Significant k-mers were found for the phenotype: " + pheno)
        print("Results are written to: " + args.out_path)
        
        if not os.path.exists(args.out_path):
            os.makedirs(args.out_path)
        kmers_pass_threshold.to_csv(args.out_path + "/" + pheno + "_kmers_list.txt", index=False, sep="\t", header=None)
    
        # =================================================================================================
        #     Convert to FASTA (based on: https://www.biostars.org/p/271977/)
        # =================================================================================================
    
        #File input
        fileInput = open(args.out_path + "/" + pheno + "_kmers_list.txt", "r")
    
        #File output
        fileOutput = open(args.out_path + "/" + pheno + "_kmers_list.fa", "w")
    
        #Seq count
        count = 1
    
        #Loop through each line in the input file
        print("Converting to FASTA...")
        for strLine in fileInput:
        
            #Strip the endline character from each input line
            strLine = strLine.rstrip("\n")
    
            #Output the header
            fileOutput.write(">" + "kmer" + str(count) + "\n")
            fileOutput.write(strLine + "\n")
    
            count = count + 1
        print ("Done.")
    
        #Close the input and output file
        fileInput.close()
        fileOutput.close()