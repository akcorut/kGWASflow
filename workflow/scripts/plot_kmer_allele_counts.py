import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats
import matplotlib.ticker as ticker
from matplotlib.pyplot import hist

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--in_file', help='', type=str, required=True)
parser.add_argument('-o', '--out_file', help='', type=str, required=True)

args = parser.parse_args()

## Read kmers_to_use.shareness file.
# This file contains the information about
# how many k-mers appeared in each individual/sample. 
# More info: https://github.com/voichek/kmersGWAS/blob/master/manual.pdf
kmers_shareness = pd.read_csv(args.in_file, 
                                sep="\t", header=None, skiprows=1)


## Plot the stats from kmers_to_use.shareness file
# x-axis is allele count and y-axis is no. of k-mers
# For every n accession/sample, the number of k-mers 
# appeared in exactly that number of accessions.
fig_dims = (10, 10)
fig, ax = plt.subplots(figsize=fig_dims)
sns_plot = sns.barplot(x=kmers_shareness[0], y=kmers_shareness[1], color="lightgreen", saturation=.5,)
sns_plot.set_xlabel("Allele count",fontsize=18, weight='bold')
sns_plot.set_ylabel("#k-mers",fontsize=18, weight='bold')
sns_plot.set_yscale("log")
# sns_plot.set(ylim=(50000, 1000000000))
sns_plot.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
sns_plot.xaxis.set_major_locator(ticker.MultipleLocator(base=5))
fig.tight_layout()
## Save the plot
sns_plot.figure.savefig(args.out_file, dpi=600)
