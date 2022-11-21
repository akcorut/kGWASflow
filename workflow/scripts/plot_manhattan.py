# import pandas and qmplot functions
import csv
import os, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from qmplot import manhattanplot
import natsort
import seaborn as sns

if __name__ == "__main__":

    ## Logging
    with open(snakemake.log[0], "w") as f:
      sys.stderr = sys.stdout = f

      ## Read the alignment data 
      align_kmers_sam = pd.read_table(snakemake.input[0], 
                                      sep='\t', comment='@', header=None,
                                      usecols=[0,2,3], names=['kmer_id', 'chr', 'bp'])

      ## Preparing the data for plotting
      align_kmers_sam['kmer'] = align_kmers_sam['kmer_id'].str.split('_').str[0]
      align_kmers_sam['p_value'] = align_kmers_sam['kmer_id'].str.split('_').str[1]
      align_kmers_sam['p_value'] = align_kmers_sam['p_value'].astype(float)
      align_kmers_sam['bp'] = align_kmers_sam['bp'].astype(int)
      
      # Sort the data by chromosome and chromosome position
      align_kmers_sam_sorted = align_kmers_sam.sort_values(by=["chr", "bp"], key=natsort.natsort_keygen())

      # Get colors for manhattan plot
      colors = sns.color_palette("colorblind").as_hex()

      # Plotting the manhattan plot
      print("Plotting...")
      f, ax = plt.subplots(figsize=(18, 9), facecolor='w', edgecolor='k')
      manhattanplot(data=align_kmers_sam_sorted,
                    snp="kmer_id",
                    chrom="chr",
                    color=colors,
                    pos="bp",
                    pv="p_value",
                    xticklabel_kws={"rotation": "vertical"},
                    ax=ax,
                    s = snakemake.params["point_size"],
                    clip_on=False)
      f.suptitle('k-mer Based GWAS Manhattan Plot for ' + snakemake.params["pheno"], fontsize=22)
      plt.xlabel('Chromosome', fontsize=18)
      plt.ylabel(r"$-log_{10}{(P)}$", fontsize=18) 
      plt.tight_layout()
      
      ## Saving the plot as pdf
      print("Plotting is done. Saving the plot...")
      plt.savefig(snakemake.output["manhattan_plot"])
