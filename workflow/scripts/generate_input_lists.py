import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--dir_path', help='', type=str, required=True)
args = parser.parse_args()

dir_names = os.listdir(args.dir_path)

for dir in dir_names:
    input_list=[]
    df=[]
    # input_list = glob.glob(args.dir_path + '/' + dir + '/' + "*.gz")
    types = ["*.gz", "*.fq", "*.fastq"]
    for type in types:
        temp = glob.glob(args.dir_path + '/' + dir + '/' + type)
        input_list += temp
    df = pd.DataFrame(sorted(input_list),columns=['input_files'])
    df.to_csv(args.dir_path + '/' + dir + '/' + "input_files.txt", index=False, header=False)

    try:
        if os.path.getsize(args.dir_path + '/' + dir + '/' + "input_files.txt") > 0:
            print("input_files.txt is successfully created for the sample: " + dir)
        else:
            print("Warning: input_files.txt is empty for the sample: "  + dir)

    except OSError as e:
        print("input_files.txt cannot be created for the sample: "  + dir)       