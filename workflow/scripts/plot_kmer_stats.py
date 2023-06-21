import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats
import matplotlib.ticker as ticker

# Logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ## Gather dir (sample names as directory) names in a list
    kmc_logs_path= snakemake.params.input_path
    dir_names= os.listdir(kmc_logs_path)

    ## Define a function to generate k-mers count stats table 
    # for each accession/sample. Accession name, no. total reads and no. of 
    # unique k-mers as columns, samples are rows.
    def generate_kmers_stats_tab(dir_path, file_name, dir_names):
        """
        Generate k-mers count stats table for each accession/sample.
        """
        # Define lists to store values
        unique_kmers=[] # list to store unique k-mers
        total_reads=[] # list to store total reads
        accessions=[] # list to store accession names
        
        # For each accession/sample gather accesssion name,
        # from dir names and no. total reads and no. of unique 
        # k-mers from kmc log files.
        for dir in dir_names: # for each accession/sample
            if os.path.isfile(dir_path + '/' + dir + '/' + file_name): # if kmc log file exists
                df= pd.read_csv(dir_path + '/' + dir + '/' + file_name, delimiter = ":", header=None) # read kmc log file
                df= df.dropna() # drop rows with NaN values
                unique_kmers.append(df.iloc[6][1]) # append unique k-mers to list
                total_reads.append(df.iloc[9][1])  # append total reads to list
                accessions.append(dir) # append accession name to list

        ## Create k-mers count stats table.
        target_df = pd.DataFrame(
            {'accessions': accessions,
             'unique_kmers': unique_kmers,
             'total_reads': total_reads
            })

        ## Convert unique k-mers and total reads values to integer.
        target_df["unique_kmers"] = pd.to_numeric(target_df["unique_kmers"]) #
        target_df["total_reads"] = pd.to_numeric(target_df["total_reads"])
        return target_df # return k-mers count stats table

    # Define functions to generete plots from k-mers count stats table (above). 
    def plot_kmers_stats_scatter(target_df, out_file, plot_name):
        """
        Plot the scatter plot with regression line.
        """
        # Define x-axis and y-axis
        X=target_df["unique_kmers"] # Unique k-mers as x-axis
        Y=target_df["total_reads"] # Total reads as y-axis
        # Plot scatter plot
        plt.figure(figsize=(10, 10)) # set figure size
        plt.scatter(X, Y, alpha=0.5, s=200) 
        slope, intercept = np.polyfit(X, Y, 1) # calculate regression line
        plt.plot(X, slope*X + intercept,color="red") # plot regression line
        plt.xlabel("No. of unique kmers", fontsize=24) # add x-axis label
        plt.ylabel("No. of total reads", fontsize=24) # add y-axis label
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.title(plot_name) # add plot title
        plt.tight_layout() # tight layout
        plt.savefig(out_file, dpi=300) # save plot
        plt.close() # close plot

    # Define a function to plot the joint plot with regression line and pearson correlation coefficient.
    def plot_kmers_stats_joint(target_df, out_file, plot_name):
        """
        Plot the joint plot with regression line and pearson correlation coefficient.
        """
        # plot joint plot
        joint_plot= sns.jointplot(x="unique_kmers", y="total_reads", 
                                  data=target_df, # data
                                  kind="reg", # regression line
                                  height=10,
                                  joint_kws={'line_kws':{'color':'red'},
                                             'scatter_kws': {'s': 150}})
        
        # pearson correlation coefficient
        r, p = stats.pearsonr(target_df["unique_kmers"], target_df["total_reads"]) # calculate pearson correlation coefficient
        phantom, = joint_plot.ax_joint.plot([], [], linestyle="", alpha=0) # add pearson correlation coefficient to the legend
        legend = joint_plot.ax_joint.legend([phantom],['r={:f}, p={:f}'.format(r,p)]) # add pearson correlation coefficient to the legend
        
        # add x-axis label, y-axis label and plot title
        joint_plot.ax_joint.set_xlabel('No. of unique kmers', fontsize=24) # add x-axis label
        joint_plot.ax_joint.set_ylabel('No. of total reads', fontsize=24) # add y-axis label
        plt.figure(figsize=(10, 10)) # set figure size
        joint_plot.ax_joint.xaxis.set_tick_params(labelsize=18)
        joint_plot.ax_joint.yaxis.set_tick_params(labelsize=18)
        plt.setp(legend.get_texts(), fontsize='14') # for legend text
        joint_plot.fig.suptitle(plot_name, y = 1) # add plot title
        joint_plot.fig.tight_layout() # tight layout
        joint_plot.savefig(out_file, dpi=300) # save plot
        plt.close() # close plot
        
    # Define a function to plot the distribution of k-mer counts for canonical and non-canonical k-mers.
    def plot_kmer_dist(target_df, out_file, plot_name):
        """
        Plot the distribution of k-mer counts for canonical and non-canonical k-mers.
        """
        # Plot the histogram of k-mer counts
        plt.figure(figsize=(10, 10)) # set figure size
        plt.hist(target_df["unique_kmers"], 
                 bins=30, 
                 color="#86bf91", 
                 alpha=0.5, 
                 edgecolor="black",
                 linewidth=1.2) #
        plt.xlabel("No. of unique kmers", fontsize=24) # add x-axis label
        plt.ylabel("No. of samples", fontsize=24) # add y-axis label
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.title(plot_name, fontsize=24) # add plot title
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) # set x-axis to scientific notation
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
        
        # Add a boxed text indicating the total number of k-mers to the top right corner of the plot
        textstr = 'Total no. of k-mers:\n{:,}'.format(target_df["unique_kmers"].sum())
        # add box around the text
        plt.gca().text(0.95, 0.95, textstr, 
                       transform=plt.gca().transAxes, 
                       fontsize=16,
                       color="black",
                       verticalalignment='top', 
                       horizontalalignment='right',
                       bbox=dict(facecolor='wheat', pad=10.0),  
                       alpha=0.5)    
                            
        plt.tight_layout() # tight layout
        plt.savefig(out_file, dpi=300) # save plot
        plt.close() # close plot    

    ## Generate stats tables for both canonized and non-canonized k-mers
    kmc_canon_stats= generate_kmers_stats_tab(dir_path=kmc_logs_path, 
                                              file_name="kmc_canon.log", 
                                              dir_names=dir_names)
    kmc_non_canon_stats= generate_kmers_stats_tab(dir_path=kmc_logs_path, 
                                                  file_name="kmc_all.log", 
                                                  dir_names=dir_names)

    ## Writing out KMC stats as a tsv table
    print("Writing out KMC stats as a tsv table...")
    kmc_canon_stats.to_csv(snakemake.output.kmc_canon_stats, 
                           index=False, sep="\t")
    kmc_non_canon_stats.to_csv(snakemake.output.kmc_all_stats, 
                               index=False, sep="\t")

    ## Plot the stats
    print("Plotting total reads vs. number of unique k-mers for canonical k-mers...")
    plot_kmers_stats_scatter(target_df=kmc_canon_stats, 
                             out_file=snakemake.output.kmc_canon_plot, 
                             plot_name="Total Reads vs. Number of Unique k-mers (Canonical)")
    print("Plotting total reads vs. number of unique k-mers for non-canonical k-mers...")
    plot_kmers_stats_scatter(target_df=kmc_non_canon_stats, 
                             out_file=snakemake.output.kmc_all_plot, 
                             plot_name="Total Reads vs. Number of Unique k-mers (Non-canonical)")
    
    print("Plotting total reads vs. number of unique k-mers for canonical k-mers (joint plot)...")
    plot_kmers_stats_joint(target_df=kmc_canon_stats, 
                           out_file=snakemake.output.kmc_canon_joint_plot, 
                           plot_name="Total Reads vs. Number of Unique k-mers (Canonical)")
    print("Plotting total reads vs. number of unique k-mers for non-canonical k-mers (joint plot)...")
    plot_kmers_stats_joint(target_df=kmc_non_canon_stats, 
                           out_file=snakemake.output.kmc_all_joint_plot, 
                           plot_name="Total Reads vs. Number of Unique k-mers (Non-canonical)")
    
    print("Plotting the distribution of k-mer counts for canonical k-mers...")
    plot_kmer_dist(target_df=kmc_canon_stats, 
                   out_file=snakemake.output.kmc_canon_kmer_dist_plot,
                   plot_name="Distribution of k-mer counts (Canonical)")
    print("Plotting the distribution of k-mer counts for non-canonical k-mers...")
    plot_kmer_dist(target_df=kmc_non_canon_stats, 
                   out_file=snakemake.output.kmc_all_kmer_dist_plot,
                   plot_name="Distribution of k-mer counts (Non-canonical)")
    
    print("Done!")