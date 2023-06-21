import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import seaborn as sns
import scipy as sp
from scipy import stats
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist

# Logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ## Read kmers_to_use.shareness file.
    # This file contains the information about
    # how many k-mers appeared in each individual/sample. 
    # More info: https://github.com/voichek/kmersGWAS/blob/master/manual.pdf
    kmers_shareness = pd.read_csv(snakemake.input.shareness, 
                                    sep="\t", header=None, skiprows=1)
    print(kmers_shareness)

    ## Plot the stats from kmers_to_use.shareness file
    # x-axis is allele count and y-axis is no. of k-mers
    # For every n accession/sample, the number of k-mers 
    # appeared in exactly that number of accessions.
    fig_dims = (10, 10)
    fig, ax = plt.subplots(figsize=fig_dims)
    sns_plot = sns.barplot(x=kmers_shareness[0], y=kmers_shareness[1], 
                           color="darksalmon", 
                           saturation=1)
    sns_plot.set_xlabel("No. of samples",fontsize=24)
    sns_plot.set_ylabel("No. of k-mers",fontsize=24)

    sns_plot.set_yscale("log")

    # set y-axis ticks
    yticks = [10**i for i in range(4, 10)]
    sns_plot.yaxis.set_major_locator(ticker.FixedLocator(yticks))

    # set y-axis tick labels
    ytick_labels = [f"$10^{i}$" for i in range(4, 10)]
    sns_plot.set_yticklabels(ytick_labels)

    # set x-axis ticks
    sns_plot.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    sns_plot.xaxis.set_major_locator(ticker.MultipleLocator(base=30))
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    # rotate x axis labels
    ax.tick_params(axis='x', rotation=45)
    fig.tight_layout()

    ## Save the plot
    sns_plot.figure.savefig(snakemake.output.kmer_allele_counts_plot, dpi=300)
