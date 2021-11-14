import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats

## Add argument parsers
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--in_dir', help='Input path for kmc log files', type=str, required=True)
parser.add_argument('-o1', '--out_table', help='Path for output tables', type=str, required=True)
parser.add_argument('-o2', '--out_plot', help='Path for output plots', type=str, required=True)
args = parser.parse_args()

## Gather dir (sample names as directory) names in a list
path= args.in_dir
dir_names= os.listdir(path)

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

## Define a function to generete plots from k-mers
# count stats table (above). This function generates
# two plots, one is a scatter plot for only one of the 
# KMC runs (canonized or non-cononized) and the other 
# plot is a joint scatter plot that uses both KMC runs.
# This funtion also fits a line.
def plot_kmers_stats(target_df, out_path, plot_name):
    # Plot the single plot
    X=target_df["unique_kmers"]
    Y=target_df["total_reads"]
    plt.figure(figsize=(9, 9))
    plt.scatter(X, Y, alpha=0.5)
    slope, intercept = np.polyfit(X, Y, 1)
    plt.plot(X, slope*X + intercept,color="red")
    plt.xlabel("No. of unique kmers")
    plt.ylabel("No. of total reads")
    plt.title(plot_name) 
    plt.savefig(out_path + "/" + plot_name + '.scatter_plot.pdf')

    # Plot the joint plot.
    sns_plot_2= sns.jointplot(x="unique_kmers", y="total_reads", 
    data=target_df, kind="reg", height=9, joint_kws={'line_kws':{'color':'red'}})
    r, p = stats.pearsonr(target_df["unique_kmers"], target_df["total_reads"])
    phantom, = sns_plot_2.ax_joint.plot([], [], linestyle="", alpha=0)
    sns_plot_2.ax_joint.legend([phantom],['r={:f}, p={:f}'.format(r,p)])
    plt.suptitle(plot_name)
    plt.savefig(out_path + "/" + plot_name + '.joint_plot.pdf')

## Generate stats tables for both canonized and non-canonized k-mers
kmc_canon_stats= generate_kmers_stats_tab(dir_path=path, file_name="kmc_canon.log", dir_names=dir_names)
kmc_non_canon_stats= generate_kmers_stats_tab(dir_path=path, file_name="kmc_all.log", dir_names=dir_names)

## Writing out KMC stats as a tsv table
kmc_canon_stats.to_csv(args.out_table + '/' + "kmc_canon.stats.tsv", index=False, sep="\t")
kmc_non_canon_stats.to_csv(args.out_table + '/' + "kmc_all.stats.tsv", index=False, sep="\t")

## Plot the stats
plot_kmers_stats(target_df=kmc_canon_stats, out_path=args.out_plot, plot_name="kmc_canon_total_reads_vs_unique_kmers")
plot_kmers_stats(target_df=kmc_non_canon_stats, out_path=args.out_plot, plot_name="kmc_all_total_reads_vs_unique_kmers")