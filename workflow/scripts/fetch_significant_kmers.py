import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import sys

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ## Set float presicion for the p-values
    pd.set_option("display.precision", 8)
    
    ## Set the threshold for k-mers GWAS results (%5 or %10 family-wise error rate)
    ## 5 is the default
    threshold = snakemake.params.threshold
    
    ## Set input and output paths
    kmers_gwas_dir = snakemake.params.in_prefix
    outdir = snakemake.output[0]
    check_file = snakemake.params.check_file_prefix + '/NO_KMERS_PASS_5PER_THRESHOLD_FOUND' 
    
    if os.path.exists(check_file):
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        with open(outdir + '/NO_KMERS_PASS_5PER_THRESHOLD_FOUND' , 'w') as fp:
            pass
        print("No k-mers passing " + str(threshold) + " percent threshold were found in any of the phenotypes.")
        print("Exiting.")
        
    else:
        ## Get the phenotype names using k-mers GWAS results directory path
        phenos=[]
        for file in os.listdir(kmers_gwas_dir):
            d = os.path.join(kmers_gwas_dir, file)
            if os.path.isdir(d):
                phenos.append(os.path.basename(d))
    
        # print(phenos) # For debugging only
    
        ## For each phenotype do the followings:
        for pheno in phenos:
        
            ## Check the threshold and choose the k-mers GWAS result 
            # based on the threshold (pass_threshold_5per or pass_threshold_10per)
            if threshold == 5:
                kmers_pass_threshold= pd.read_csv(kmers_gwas_dir + '/' + pheno + '/kmers/pass_threshold_5per',
                                                  sep="\t")
            else:
                kmers_pass_threshold= pd.read_csv(kmers_gwas_dir + '/' + pheno + '/kmers/pass_threshold_10per',
                                                  sep="\t")
    
            ## Check if any significant k-mers associated with given phenotype were found
            if kmers_pass_threshold.empty:
                pass
                print("No k-mers passing " + str(threshold) + " percent threshold were found for the phenotype: " + pheno)
                print("")
    
            else:
                print("Signifincat k-mers passing " + str(threshold)  + " percent threshold were found for the phenotype: " + pheno)
    
                ## Get significant k-mers without k-mer number
                kmers_pass_threshold['rs']= kmers_pass_threshold['rs'].apply(lambda x: x.split('_')[0])
    
                # print(kmers_pass_threshold) # For debugging only
    
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
    
                ## Write significant k-mers list in TEXT format 
                print("Writing significant k-mers list in TEXT format...")
                kmers_pass_threshold['rs'].to_csv(outdir + "/" + pheno + "_kmers_list.txt", 
                                                 index=False, sep="\t", header=None)
                print ("Done.")
    
                ## Make a new column with k-mer number and its p_value (e.g. kmer1_8.543334e-11)
                kmers_pass_threshold['kmer_name']= 'kmer' + (kmers_pass_threshold.index+1).astype(str) + '_' + kmers_pass_threshold['p_lrt'].astype(str)
    
                # print(kmers_pass_threshold) # For debugging only
    
                ### Write significant k-mers list in FASTA format (based on: https://www.biostars.org/p/271977/)
                ## Fasta file output
                fileOutput = open(outdir + "/" + pheno + "_kmers_list.fa", "w")
    
                ## Seq count
                count = 1
    
                ## Loop through each line in the input file
                print("Writing significant k-mers list in FASTA format...")
                for index, row in kmers_pass_threshold.iterrows():
                
                    #Output the header
                    fileOutput.write(">" + row['kmer_name'] + "\n")
                    fileOutput.write(row['rs'] + "\n")
    
                    count = count + 1
                print ("Done.")
                print("")
    
                #Close the input and output file
                fileOutput.close()