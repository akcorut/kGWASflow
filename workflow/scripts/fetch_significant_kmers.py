import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import sys

# Logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ## Set float presicion for the p-values
    pd.set_option("display.precision", 8)
    
    # Input directory path for kmersGWAS results
    res_dir_path = os.path.dirname(snakemake.params.in_prefix)

    ## Set the threshold for k-mers GWAS results (%5 or %10 family-wise error rate)
    ## 5 is the default
    threshold = snakemake.params.threshold
    
    # If threshold is not 5 or 10, print an error message and exit
    if threshold != 5 and threshold != 10:
        print("ERROR: threshold value must be 5 or 10. \nPlease enter a valid threshold (5 or 10). \nExiting...")
        sys.exit(1)
    
    # Define the output directory
    outdir = snakemake.output.dir
    print("Output directory: " + outdir)
    
    # Define a function to fetch significant k-mers
    def fetch_kmers(input_dir, pheno, out_dir, threshold):
        """
        Fetch significant k-mers for a given phenotype from the k-mers GWAS results.
        """

        kmers_pass_threshold = pd.read_csv(input_dir + '/' + pheno +  '/kmers/pass_threshold_{threshold}per'.
                                           format(threshold=threshold),
                               sep="\t") # read the k-mers that pass the threshold
        
        if not os.path.exists(out_dir):
                os.makedirs(out_dir) # create the output directory if it does not exist
                
        if kmers_pass_threshold.empty: # if no k-mers pass the 5% threshold:
            # create a file to indicate that no k-mers pass the 5% threshold
            with open(out_dir + '/{pheno}_NO_KMERS_PASS_5PER_THRESHOLD_FOUND'.
                      format(pheno=pheno), 'a'): 
                print("No k-mers pass the {threshold}% threshold for the phenotype {pheno}.".
                      format(threshold=threshold, pheno=pheno)) # print a message to the user and exit the script
                pass
            
        else: # if k-mers pass thethreshold
            
        ## Get significant k-mers without k-mer number
            kmers_pass_threshold['rs']= kmers_pass_threshold['rs'].apply(lambda x: x.split('_')[0])
    
            # print(kmers_pass_threshold) # For debugging only
    
            ## Write significant k-mers list in TEXT format 
            print("Writing significant k-mers list in TEXT format...")
            kmers_pass_threshold['rs'].to_csv(out_dir + "/" + pheno + "_kmers_list.txt", 
                                             index=False, 
                                             sep="\t", 
                                             header=None) # write the k-mers list in a text file
            print ("Done.")
    
            ## Make a new column with k-mer number and its p_value (e.g. kmer1_8.543334e-11)
            kmers_pass_threshold['kmer_name']= 'kmer' + (kmers_pass_threshold.index+1).astype(str) + '_' + kmers_pass_threshold['p_lrt'].astype(str)
    
            ### Write significant k-mers list in FASTA format (based on: https://www.biostars.org/p/271977/)
            ## Fasta file output
            fileOutput = open(out_dir + "/" + pheno + "_kmers_list.fa", "w")
    
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

    # Get the list of phenotypes
    phenos = [os.path.basename(x) for x in glob.glob(res_dir_path + '/*')]
    print("Phenotypes: " + str(phenos)) # print the list of phenotypes
    
    for pheno in phenos: # for each phenotype
        
        print("Fetching significant k-mers for phenotype: " + pheno)
        print("Threshold: " + str(threshold) + "%")
        
        # Fetch significant k-mers for the current phenotype
        fetch_kmers(input_dir=res_dir_path, pheno=pheno, out_dir=outdir , threshold=threshold)
        
    print("Done!.")


 
