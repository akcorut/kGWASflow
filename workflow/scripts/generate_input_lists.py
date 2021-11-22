import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import sys

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    dir_path = snakemake.params.prefix # Directory path of input reads
    
    # Get dir names (each dir name shoud be same as sample names)
    dir_names = os.listdir(dir_path)

    # For each dir/sample create a list 
    # that contains full path of each 
    # fastq/fq/.gz files exist in that folder
    for dir in dir_names:
        input_list=[]
        df=[]
        types = ["*.gz", "*.fq", "*.fastq"]
        for type in types:
            temp = glob.glob(dir_path + '/' + dir + '/' + type)
            input_list += temp
        df = pd.DataFrame(sorted(input_list),columns=['input_files'])

        # Write out list of fastq files to input_files.txt
        # and save it in current dir/sample
        df.to_csv(dir_path + '/' + dir + '/' + "input_files.txt", index=False, header=False)

        # Check if input_files.txt is created and not empty
        try:
            if os.path.getsize(dir_path + '/' + dir + '/' + "input_files.txt") > 0:
                print("input_files.txt is successfully created for the sample: " + dir)
            else:
                print("Warning: input_files.txt is empty for the sample: "  + dir)

        except OSError as e:
            print("input_files.txt cannot be created for the sample: "  + dir)       