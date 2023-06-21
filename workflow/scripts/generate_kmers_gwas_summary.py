import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

# Logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
    # Read the kmersGWAS results file (pass_threshold_5per)
    pass_threshold_5per = pd.read_csv(snakemake.input.kmers_gwas_res,
                                      sep="\t") 

    # Input directory path for kmersGWAS results
    res_dir_path = os.path.dirname(snakemake.input.kmers_gwas_res)

    # Read the pass_threshold_10per file
    pass_threshold_10per = pd.read_csv(res_dir_path + '/pass_threshold_10per',
                                       sep="\t")

    # Get the permutation based p-value thresholds values for 5% and 10% 
    threshold_5per = pd.read_csv(res_dir_path + '/threshold_5per',
                                 sep="\t", 
                                 names=['threshold'])
    threshold_10per = pd.read_csv(res_dir_path + '/threshold_10per', 
                                 header=None, 
                                 sep="\t", 
                                 names=['threshold'])

    # read the phenotype_value.assoc.txt.gz file
    pheno_val_assoc = pd.read_csv(res_dir_path + '/output/phenotype_value.assoc.txt.gz', 
                              compression='gzip',
                              sep="\t")

    # Phenotype name
    pheno= snakemake.params.pheno

    # Number of k-mers filtered
    kmers_number=snakemake.params.kmers_number

    # Define a function to generate a table with the k-mers, p-values and p-value thresholds
    def generate_kmers_gwas_results_table(pass_threshold, threshold, pheno, output_file):
        
        # Generate an empty dataframe with the kmer, p-value and threshold
        res_table = pd.DataFrame(columns=['kmers_pass_threshold', 'pval', 'pval_threshold'])
        
        # Fill the dataframe with the kmer, p-value, threshold, phenotype, -log10(p-value) and allele frequency
        res_table['kmers_pass_threshold'] = pass_threshold['rs'].apply(lambda x: x.split('_')[0]) # get the kmer from the rs column

        # Get the p-value from the p_lrt column and round it to 6 decimal places
        res_table['pval'] = pass_threshold['p_lrt']
        
        # Get the -log10(p-value) from the p-value column
        res_table['-log10(pval)'] = -np.log10(res_table['pval'])
        
        res_table['pval'] = res_table['pval'].map('{:.6e}'.format) # round to 6 decimal places
        res_table['-log10(pval)'] = res_table['-log10(pval)'].map('{:.5f}'.format) # round to 5 decimal places

        # Get the allele frequency from the af column 
        res_table['allele_freq'] = pass_threshold['af'] # add the allele frequency column

        # Get the p-value threshold
        res_table['pval_threshold'] = threshold['threshold'][0] # add the p-value threshold column
        res_table['phenotype'] = pheno # add the phenotype column

        print(res_table.head())
        # Save the table to a tsv file
        print('Saving the results summary table to a tsv file...')
        res_table.to_csv(output_file, index=False, 
                     sep='\t', 
                     columns=['kmers_pass_threshold', 'phenotype', 'allele_freq', 'pval', '-log10(pval)', 'pval_threshold'])    
        print('Generating the results summary table is done!')
        
    # Define a function to plot the histogram of the -log10(p-value) values
    def plot_pval_histogram(pheno_val_assoc, threshold_5per, threshold_10per, kmers_number, pheno, out_dir): 
        
        # If output directory does not exist, create it
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Plot the histogram of the -log10(p-value) values
        f, ax = plt.subplots(figsize=(16, 8), facecolor='w', edgecolor='k')
        plt.hist(-np.log10(pheno_val_assoc['p_lrt']), bins=50, color='g', alpha=0.5)
        
        # Add labels to the axes
        plt.xlabel('-log10(pval)', fontsize=24)
        plt.ylabel('Frequency', fontsize=24)

        # Add a vertical line to the histogram based on the 5% threshold
        plt.axvline(x=threshold_5per['threshold'][0], color='r', linestyle='--')
         # Add a caption to the right bottom corner of the outside of the plot
        ax.text(1.0, -0.1, "*Red dashed line indicates 5% family-wise error-rate threshold", color='r', transform=ax.transAxes, ha='right', va='bottom', fontsize=12)
        
        # Add label to the top of the vertical line
        ax.text(threshold_5per['threshold'][0]+0.3, 1, '5% threshold', 
                rotation=0, color='r', va='top' , ha='left', 
                bbox=dict(facecolor='w', edgecolor='r', boxstyle='round,pad=0.5'),
                transform=ax.get_xaxis_transform(),
                fontsize=12)
        
        # Add a vertical line to the histogram based on the 10% threshold
        plt.axvline(x=threshold_10per['threshold'][0], color='b', linestyle='--')
        
        # Add a caption to the right bottom corner of the outside of the plot
        ax.text(1.0, -0.13, "*Blue dashed line indicates 10% family-wise error-rate threshold", color='b', transform=ax.transAxes, ha='right', va='bottom', fontsize=12)
        # Add label to the top of the vertical line
        ax.text(threshold_10per['threshold'][0]-0.3, 1, '10% threshold', 
                rotation=0, color='b', va='top' , ha='right', 
                bbox=dict(facecolor='w', edgecolor='b', boxstyle='round,pad=0.5'),
                transform=ax.get_xaxis_transform(),
                fontsize=12)

        # Add title to the histogram
        plt.title('Histogram of (-log10) p-values - k-mers passing the first filtering step (n={kmers_num})'.format(kmers_num=kmers_number), 
                  fontsize=16, y=1.02)
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        
        # Adjust the layout of the plot
        plt.tight_layout() 
        print('Saving the histogram plot to a pdf file...')
        
        # Save the plot to a PDF file
        plt.savefig(out_dir + "/{pheno}.kmers_pass_first_filter.pval_hist.pdf".format(pheno=pheno))
        print("Generating the p-value histogram plot is done!")
    
    outdir = os.path.dirname(snakemake.output.res_sum_5per)
    
    # Check if there are k-mers that pass the 5% threshold
    if pass_threshold_5per.empty: # if no k-mers pass the 5% threshold
        with open(outdir + "/NO_KMERS_PASS_5PER_THRESHOLD_FOUND", 'a'): # create a file to indicate that no k-mers pass the 5% threshold
            print("No k-mers pass the 5% threshold.") # print a message to the user

    # generate the table with the k-mers that pass the 5% threshold
    print("Generating kmersGWAS summary table with k-mers that pass the 5% threshold...")
    generate_kmers_gwas_results_table(pass_threshold=pass_threshold_5per, threshold=threshold_5per, pheno=pheno, output_file=snakemake.output.res_sum_5per) 

    # Check if there are k-mers that pass the 10% threshold
    if pass_threshold_10per.empty: # if no k-mers pass the 10% threshold
        with open(outdir + "/NO_KMERS_PASS_10PER_THRESHOLD_FOUND", 'a'): # create a file to indicate that no k-mers pass the 10% threshold
            print("No k-mers pass the 10% threshold.") # print a message to the user and exit the script

    # generate the table with the k-mers that pass the 10% threshold
    print("Generating kmersGWAS summary table with k-mers that pass the 10% threshold...")
    generate_kmers_gwas_results_table(pass_threshold=pass_threshold_10per, 
                                      threshold=threshold_10per, 
                                      pheno=pheno, 
                                      output_file=snakemake.output.res_sum_10per) 
    
    # if k-mers pass the 5% threshold or the 10% threshold 
    if not pass_threshold_5per.empty or not pass_threshold_10per.empty:
    # Plot the histogram of the -log10(p-value) values
        print("Plotting the histogram of the -log10(p-value) values...")
        plot_pval_histogram(pheno_val_assoc=pheno_val_assoc, 
                            threshold_5per=threshold_5per, 
                            threshold_10per=threshold_10per,
                            kmers_number=kmers_number, 
                            pheno=pheno,
                            out_dir=snakemake.output.pval_hist_dir)
    # if no k-mers pass the 5% threshold and no k-mers pass the 10% threshold
    if pass_threshold_5per.empty and pass_threshold_10per.empty:
        # If output directory does not exist, create it
        if not os.path.exists(snakemake.output.pval_hist_dir):
            os.makedirs(snakemake.output.pval_hist_dir)
        
        # Create a file to indicate that no k-mers pass the 5% threshold
        with open(snakemake.output.pval_hist_dir + "/NO_KMERS_PASS_10PER_THRESHOLD_FOUND", 'a'):
            print("No k-mers pass the 5% or the 10% threshold. No plots will be generated.") # print a message to the user and exit the script
    
    print("Done!")