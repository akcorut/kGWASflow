import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ## Gather dir (sample names as directory) names in a list
    kmc_logs_path= snakemake.params.input_path
    dir_names= os.listdir(kmc_logs_path)

    ## Define a function to generate k-mers count stats table 
    # for each accession/sample. Accession name, no. total reads and no. of 
    # unique k-mers as columns, samples are rows.
    def generate_kmers_stats_tab(dir_path, file_name, dir_names):
        unique_kmers=[]
        total_reads=[]
        accessions=[]
        ## For each accession/sample gather accesssion name,
        # from dir names and no. total reads and no. of unique 
        # k-mers from kmc log files.
        for dir in dir_names:
            if os.path.isfile(dir_path + '/' + dir + '/' + file_name):
                df= pd.read_csv(dir_path + '/' + dir + '/' + file_name, delimiter = ":", header=None)
                df= df.dropna()
                unique_kmers.append(df.iloc[6][1])
                total_reads.append(df.iloc[9][1])
                accessions.append(dir)

        ## Create k-mers count stats table.
        target_df = pd.DataFrame(
            {'accessions': accessions,
             'unique_kmers': unique_kmers,
             'total_reads': total_reads
            })

        ## Convert unique k-mers and total reads values to integer.
        target_df["unique_kmers"] = pd.to_numeric(target_df["unique_kmers"])
        target_df["total_reads"] = pd.to_numeric(target_df["total_reads"])
        return target_df

    ## Define functions to generete plots from k-mers
    ## count stats table (above). 
    def plot_kmers_stats_scatter(target_df, out_file, plot_name):
        # Plot the single plot
        X=target_df["unique_kmers"]
        Y=target_df["total_reads"]
        plt.figure(figsize=(10, 10))
        plt.scatter(X, Y, alpha=0.5)
        slope, intercept = np.polyfit(X, Y, 1)
        plt.plot(X, slope*X + intercept,color="red")
        plt.xlabel("No. of unique kmers")
        plt.ylabel("No. of total reads")
        plt.title(plot_name) 
        plt.savefig(out_file)

    def plot_kmers_stats_joint(target_df, out_file, plot_name):
        # Plot the joint plot.
        sns_plot_2= sns.jointplot(x="unique_kmers", y="total_reads", 
            data=target_df, kind="reg", height=10, 
            joint_kws={'line_kws':{'color':'red'}})
        r, p = stats.pearsonr(target_df["unique_kmers"], target_df["total_reads"])
        phantom, = sns_plot_2.ax_joint.plot([], [], linestyle="", alpha=0)
        sns_plot_2.ax_joint.legend([phantom],['r={:f}, p={:f}'.format(r,p)])
        sns_plot_2.ax_joint.set_xlabel('No. of unique kmers')
        sns_plot_2.ax_joint.set_ylabel('No. of total reads')
        sns_plot_2.fig.suptitle(plot_name, y = 1)
        sns_plot_2.fig.tight_layout() 
        plt.savefig(out_file)

    ## Generate stats tables for both canonized and non-canonized k-mers
    kmc_canon_stats= generate_kmers_stats_tab(dir_path=kmc_logs_path, file_name="kmc_canon.log", dir_names=dir_names)
    kmc_non_canon_stats= generate_kmers_stats_tab(dir_path=kmc_logs_path, file_name="kmc_all.log", dir_names=dir_names)

    ## Writing out KMC stats as a tsv table
    kmc_canon_stats.to_csv(snakemake.output.kmc_canon_stats, index=False, sep="\t")
    kmc_non_canon_stats.to_csv(snakemake.output.kmc_all_stats, index=False, sep="\t")

    ## Plot the stats
    plot_kmers_stats_scatter(target_df=kmc_canon_stats, out_file=snakemake.output.kmc_canon_plot, plot_name="Total Reads vs. Number of Unique k-mers (Canonical)")
    plot_kmers_stats_scatter(target_df=kmc_non_canon_stats, out_file=snakemake.output.kmc_all_plot, plot_name="Total Reads vs. Number of Unique k-mers (Non-canonical)")
    
    plot_kmers_stats_joint(target_df=kmc_canon_stats, out_file=snakemake.output.kmc_canon_joint_plot, plot_name="Total Reads vs. Number of Unique k-mers (Canonical)")
    plot_kmers_stats_joint(target_df=kmc_non_canon_stats, out_file=snakemake.output.kmc_all_joint_plot, plot_name="Total Reads vs. Number of Unique k-mers (Non-canonical)")