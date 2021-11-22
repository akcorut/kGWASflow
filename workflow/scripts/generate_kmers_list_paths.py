import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import sys

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    in_dir = snakemake.params.input_dir # Path for k-mers counts folder
    samp_tab = snakemake.input.samples_tab # Path for sampels sheet (samples.tsv)
    
    dir_names = os.listdir(in_dir)

    ## Read samples.tsv (provided by the user) file 
    # to get the order of samples
    samples_tab = pd.read_csv(samp_tab, 
                                    sep="\t",
                                    usecols=["sample_name"],
                                    dtype='object')

    samples_tab = samples_tab.drop_duplicates()
    ## Loop through kmers_count results for each
    # accession/sample directory and make a list of 
    # kmers_with_strand file paths for each accession/sample
    input_list=[]
    dirs=[]
    for dir in dir_names:
        if os.path.isfile(in_dir + '/' + dir + '/' + "kmers_with_strand"):
            input_list.append(in_dir + '/' + dir + '/' + "kmers_with_strand")
            dirs.append(dir)

    ## Make a dataframe with columns accession name and 
    # kmers_with_strand file path for that acession
    list_paths = pd.DataFrame(input_list,columns=['input_files'])
    list_paths["accession"] = dirs

    ## Make sure the output file has same sample order as in samples.tsv file
    list_paths_sorted = samples_tab.merge(list_paths, left_on='sample_name',
              right_on = 'accession').drop('sample_name',axis=1)

    ## Write out output file
    list_paths_sorted.to_csv(snakemake.params.out_dir + '/' + "kmers_list_paths.txt", index=False, header=False, sep="\t")
