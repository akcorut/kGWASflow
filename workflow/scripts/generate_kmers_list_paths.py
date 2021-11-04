import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--in_dir', help='', type=str, required=True)
parser.add_argument('-o', '--out_dir', help='', type=str, required=True)
args = parser.parse_args()

dir_names = os.listdir(args.in_dir)

input_list=[]
dirs=[]
for dir in dir_names:
    if os.path.isfile(args.in_dir + '/' + dir + '/' + "kmers_with_strand"):
        input_list.append(args.in_dir + '/' + dir + '/' + "kmers_with_strand")
        dirs.append(dir)
print(dirs)

df = pd.DataFrame(input_list,columns=['input_files'])
df["accession"] = dirs
df.to_csv(args.out_dir + '/' + "kmers_list_paths.txt", index=False, header=False, sep="\t")
print(input_list)