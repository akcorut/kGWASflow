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

      # Make a column of minus log10 p-values
      align_kmers_sam_sorted['minuslog10pvalue'] = -np.log10(align_kmers_sam_sorted.p_value)
      
      ## Get min & max minus log10 p-values for y axis limits
      y_max = align_kmers_sam_sorted['minuslog10pvalue'].max()
      y_min = align_kmers_sam_sorted['minuslog10pvalue'].min()

      ## Check if only one chromosome is provided for the manhattan plot
      num_of_chrs = len(pd.unique(align_kmers_sam_sorted['chr']))
      
      # Plotting the manhattan plot
      print("Plotting...")
      
      ## If only one chromosome is provided, plot the k-mer's position on
      ## that chromosome on the x axis
      if num_of_chrs == 1:
        f, ax = plt.subplots(figsize=(18, 9), facecolor='w', edgecolor='k')
        manhattanplot(data=align_kmers_sam_sorted,
                      snp="kmer_id",
                      chrom="chr",
                      CHR= pd.unique(align_kmers_sam_sorted['chr']),
                      pos="bp",
                      pv="p_value",
                      suggestiveline=None,  # Turn off suggestiveline
                      genomewideline=None,  # Turn off genomewideline
                      xticklabel_kws={"rotation": "vertical"},
                      ax=ax,
                      s = snakemake.params["point_size"],
                      clip_on=False)
        ax.set_ylim([y_min-5, y_max+5]) # Set y axis limits
        f.suptitle('k-mer Based GWAS Manhattan Plot for ' + snakemake.params["pheno"], fontsize=22)
        plt.xlabel('Chromosome: ' + pd.unique(align_kmers_sam_sorted['chr'])[0], fontsize=18)
        plt.ylabel(r"$-log_{10}{(P)}$", fontsize=18) 
        plt.tight_layout()
      
      ## If more than one chromosome is provided, use all chromosomes
      if num_of_chrs > 1:
        f, ax = plt.subplots(figsize=(18, 9), facecolor='w', edgecolor='k')
        manhattanplot(data=align_kmers_sam_sorted,
                      snp="kmer_id",
                      chrom="chr",
                      color=colors,
                      pos="bp",
                      pv="p_value",
                      suggestiveline=None,  # Turn off suggestiveline
                      genomewideline=None,  # Turn off genomewideline
                      xticklabel_kws={"rotation": "vertical"},
                      ax=ax,
                      s = snakemake.params["point_size"],
                      clip_on=False)
        ax.set_ylim([y_min-5, y_max+5]) # Set y axis limits
        f.suptitle('k-mer Based GWAS Manhattan Plot for ' + snakemake.params["pheno"], fontsize=22)
        plt.xlabel('Chromosome', fontsize=18)
        plt.ylabel(r"$-log_{10}{(P)}$", fontsize=18) 
        plt.tight_layout()
      
      ## Saving the plot as pdf
      print("Plotting is done. Saving the plot...")
      plt.savefig(snakemake.output["manhattan_plot"])
