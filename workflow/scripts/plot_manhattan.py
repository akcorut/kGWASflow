# import pandas and qmplot functions
import csv
import os, glob, shutil
import re
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from qmplot import manhattanplot
import natsort
from natsort import natsorted
import seaborn as sns
import pysam

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
      
      # Drop na values
      align_kmers_sam = align_kmers_sam.dropna(subset=['p_value'])

      # Sort the data by chromosome and chromosome position
      align_kmers_sam_sorted = align_kmers_sam.sort_values(by=["chr", "bp"], key=natsort.natsort_keygen())

      # Get colors for manhattan plot
      colors = sns.color_palette("colorblind").as_hex()
      colors_2 = sns.color_palette("husl").as_hex()

      # Make a column of minus log10 p-values
      align_kmers_sam_sorted['minuslog10pvalue'] = -np.log10(align_kmers_sam_sorted.p_value)
      
      ## Get min & max minus log10 p-values for y axis limits
      y_max = align_kmers_sam_sorted['minuslog10pvalue'].max()
      y_min = align_kmers_sam_sorted['minuslog10pvalue'].min()
      print("y_max: " + str(y_max))
      print("y_min: " + str(y_min))
      ## Check if only one chromosome is provided for the manhattan plot
      num_of_chrs = len(pd.unique(align_kmers_sam_sorted['chr']))
      
      # Define a function to extract chromosome names and lengths from a SAM file header
      def extract_chromosome_info(sam_file):
        """
        Extract chromosome names and lengths from a SAM file header and return as a Pandas DataFrame with columns "chr" and "bp".
        """
        chromosome_info = {} # dictionary to store chromosome names and lengths
        pattern = re.compile(r'^([Cc][Hh][Rr])?\d*[XYM]?$') # chromosome names pattern
        
        with pysam.AlignmentFile(sam_file, "r") as sam: # open SAM file
            for header_line in sam.header["SQ"]: # iterate over SQ header lines
                chromosome_name = header_line["SN"] # get chromosome name
                length = header_line["LN"] # get chromosome length
                if chromosome_name.startswith("chr"): # if chromosome name starts with "chr"
                    name = chromosome_name # use name as is
                else: # if chromosome name does not start with "chr"
                    match = pattern.match(chromosome_name) # 
                    if match: # if chromosome name matches pattern
                        # Convert name to "chrX" format
                        name = "chr" + match.group(1)
                    else: # if chromosome name does not match pattern
                        continue # skip chromosome
                chromosome_info[name] = length # add chromosome name and length to dictionary
                    
        # Convert dictionary to DataFrame
        df = pd.DataFrame(chromosome_info.items(), columns=["chr", "bp"])
        
        return df # return the DataFrame
      
      # Plotting the manhattan plot
      print("Plotting...")
      
      # Set font sizes
      tick_fontsize = snakemake.params["tick_fontsize"]
      label_fontsize = snakemake.params["label_fontsize"]
      title_fontsize = snakemake.params["title_fontsize"]
      
      # Set the figure dpi
      dpi = snakemake.params["dpi"]
      figure_width = snakemake.params["figure_width"]
      figure_height = snakemake.params["figure_height"]
      
      ## If only one chromosome is provided, plot the k-mer's position on
      ## that chromosome on the x axis
      if num_of_chrs == 1:
        f, ax = plt.subplots(figsize=(figure_width, figure_height), facecolor='w', edgecolor='k')
        manhattanplot(data=align_kmers_sam_sorted,
                      snp="kmer_id",
                      chrom="chr",
                      CHR= pd.unique(align_kmers_sam_sorted['chr']),
                      color=colors_2,
                      pos="bp",
                      pv="p_value",
                      suggestiveline=None,  # Turn off suggestiveline
                      genomewideline=None,  # Turn off genomewideline
                      xticklabel_kws={"rotation": "vertical"},
                      ax=ax,
                      s = snakemake.params["point_size"],
                      clip_on=False)
        ax.set_ylim([y_min-0.5, y_max+1]) # Set y axis limits
        
        # Set x axis tick interval
        xtick_interval = snakemake.params["xtick_interval"]
        
        # Calculate the minimum and maximum of your data, rounded to the nearest multiple of xtick_interval
        min_val = align_kmers_sam_sorted['bp'].min() // xtick_interval * xtick_interval
        max_val = (align_kmers_sam_sorted['bp'].max() // xtick_interval + 1) * xtick_interval

        # Determine the range and midpoint of your data
        range_val = align_kmers_sam_sorted['bp'].max() - align_kmers_sam_sorted['bp'].min()
        midpoint = align_kmers_sam_sorted['bp'].min() + (range_val / 2)

        # Calculate the left and right limits for your x-axis
        half_width = (max_val - min_val) / 2
        left_limit = max(0, midpoint - half_width)  # ensure left_limit doesn't go below 0
        right_limit = midpoint + half_width

        # Set the x limits of your plot
        ax.set_xlim([left_limit, right_limit])

        # Calculate the closest multiple of xtick_interval that is greater than or equal to left_limit
        starting_tick = np.ceil(left_limit / xtick_interval) * xtick_interval

        # Generate the tick locations based on the updated limits and interval
        xtick_locs = np.arange(starting_tick, right_limit + 100, xtick_interval)
        
        f.suptitle('k-mer Based GWAS Manhattan Plot for ' + snakemake.params["pheno"], fontsize=title_fontsize)
        plt.xlabel('Chromosome: ' + pd.unique(align_kmers_sam_sorted['chr'])[0], fontsize=label_fontsize)
        plt.ylabel(r"$-log_{10}{(P)}$", fontsize=label_fontsize) 
        plt.xticks(xtick_locs, fontsize = tick_fontsize)
        plt.yticks(fontsize = tick_fontsize)
        plt.tight_layout()
      
      ## If more than one chromosome is provided, use all chromosomes
      if num_of_chrs > 1:

        # Extract chromosome names and lengths from the SAM file header
        chrom_info_tab = extract_chromosome_info(snakemake.input[0])
        chrom_names = natsorted(chrom_info_tab['chr'].tolist())
        
        # Add extra chromosome names and lengths to the data frame
        align_kmers_sam_with_all_chrom = pd.concat([align_kmers_sam, chrom_info_tab], ignore_index=True)
        align_kmers_sam_with_all_chrom = align_kmers_sam_with_all_chrom[align_kmers_sam_with_all_chrom['chr'].isin(chrom_names)]

        # Sort the data by chromosome and chromosome position
        align_kmers_sam_with_all_chrom_sorted = align_kmers_sam_with_all_chrom.sort_values(by=["chr", "bp"], key=natsort.natsort_keygen()) 
        # Fill NaN values with 1
        align_kmers_sam_with_all_chrom_sorted['p_value'] = align_kmers_sam_with_all_chrom_sorted['p_value'].fillna(1) 

        # Plot the dots of chromosome length rows wiht the lowest opacity 
        extra_rows = align_kmers_sam_with_all_chrom_sorted[align_kmers_sam_with_all_chrom_sorted['p_value'] == 1]

        f, ax = plt.subplots(figsize=(figure_width, figure_height), facecolor='w', edgecolor='k')
        manhattanplot(data=align_kmers_sam_with_all_chrom_sorted,
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
                      clip_on=True)
        ax.set_ylim([y_min-2, y_max+1]) # Set y axis limits
        plt.scatter(extra_rows['bp'], -np.log10(extra_rows['p_value']), alpha=0)
        
        # Calculate the cumulative distances from the start of each chromosome and store them in a list
        chrom_ends = chrom_info_tab['bp'].cumsum().tolist()

        # Plot the vertical lines for the end of each chromosome
        for end_position in chrom_ends:
          plt.axvline(x=end_position, color='grey', linestyle='--', alpha=0.2)
          
        # Add a caption to the right bottom corner of the outside of the plot
        caption = '*Vertical dashed lines indicate chromosome boundaries.'
        ax.text(1.0, -0.2, caption, transform=ax.transAxes, ha='right', va='bottom', fontsize=12)
        
        # Set the title of the plot
        f.suptitle('k-mer Based GWAS Manhattan Plot for ' + snakemake.params["pheno"], fontsize=title_fontsize) 
        
        # Set the x and y axis labels
        plt.xlabel('Chromosome', fontsize=label_fontsize) # Set the x axis label
        plt.ylabel(r"$-log_{10}{(P)}$", fontsize=label_fontsize) # Set the y axis label
        plt.xticks(fontsize = tick_fontsize)
        plt.yticks(fontsize = tick_fontsize)
        
        # Adjust the plot layout
        plt.tight_layout() 
      
      ## Saving the plot as pdf
      print("Plotting is done. Saving the plot...")
      plt.savefig(snakemake.output["manhattan_plot"], dpi=dpi)
