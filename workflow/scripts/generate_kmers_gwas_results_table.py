import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description='')
parser.add_argument('-dir', '--dir_path', help='Input directory path for kmersGWAS results', type=str, required=True)
parser.add_argument('-out', '--out_path', help='Output directory path to store results tables', type=str, required=True)
parser.add_argument('--use_10_per', default=False, action="store_true", help="Activate to get 10per_threshold results")
args = parser.parse_args()


# Define a function to parse kmersGWAS results (pass_threshold_5per or pass_threshold_10per)
def parse_kmersGWAS_results(input_dir="", threshold=5):

    # Get directory/phenotype names
    dir_names = glob.glob(input_dir + '/' + r'**', recursive=True)

    if threshold == 5:

        sig_kmers=[]
        pheno=[]
        num_sig_kmers=[]
        pheno_kmer_dict=defaultdict(list)
        # For each pass_threshold_5per file get the phenotype 
        # name and the significantly associated k-mers
        for dir in dir_names:
            if os.path.isfile(dir + "/kmers/pass_threshold_5per"):
                pass_threshold_5per = pd.read_csv(dir + '/kmers/pass_threshold_5per',sep="\t")
                pheno.append(dir.split('/')[-1])

                sig_kmers=[]
                sig_kmers= pass_threshold_5per['rs'].apply(lambda x: x.split('_')[0]).tolist()

                # Get number of significant k-mers associated with the phenotype
                num_sig_kmers.append(pass_threshold_5per.shape[0])

                # Create a phenotype - significant k-mers dictionary
                pheno_kmer_dict[dir.split('/')[-1]].append(sig_kmers)


        df_5per = pd.DataFrame(pd.Series(pheno_kmer_dict).reset_index()).set_axis(['pheno','kmers_pass_threshold_5per'],1,inplace=False)
        df_5per['kmers_pass_threshold_5per'] = df_5per['kmers_pass_threshold_5per'].apply(lambda x: x[0])
        df_5per['num_kmers_pass_threshold_5per'] = num_sig_kmers
        df_5per['num_kmers_pass_threshold_5per'] = df_5per['num_kmers_pass_threshold_5per'].apply(pd.to_numeric)
        print("{kmers_num} out of {size} phenotypes has k-mers that pass the 5per threshold.".format(kmers_num=(df_5per['num_kmers_pass_threshold_5per'] != 0).sum(), size= df_5per.shape[0]))
        return(df_5per)

    if threshold== 10:

        sig_kmers=[]
        pheno=[]
        num_sig_kmers=[]
        pheno_kmer_dict=defaultdict(list)
        for dir in dir_names:
            if os.path.isfile(dir + "/kmers/pass_threshold_10per"):
                pass_threshold_10per = pd.read_csv(dir + '/kmers/pass_threshold_10per',sep="\t")
                pheno.append(dir.split('/')[-1])

                sig_kmers=[]
                sig_kmers= pass_threshold_10per['rs'].apply(lambda x: x.split('_')[0]).tolist()

                num_sig_kmers.append(pass_threshold_10per.shape[0])

                pheno_kmer_dict[dir.split('/')[-1]].append(sig_kmers)

        df_10per = pd.DataFrame(pd.Series(pheno_kmer_dict).reset_index()).set_axis(['pheno','kmers_pass_threshold_10per'],1,inplace=False)
        df_10per['kmers_pass_threshold_10per'] = df_10per['kmers_pass_threshold_10per'].apply(lambda x: x[0])
        df_10per['num_kmers_pass_threshold_10per'] = num_sig_kmers
        df_10per['num_kmers_pass_threshold_10per'] = df_10per['num_kmers_pass_threshold_10per'].apply(pd.to_numeric)
        print("{kmers_num} out of {size} phenotypes has k-mers that pass the 10per threshold.".format(kmers_num=(df_10per['num_kmers_pass_threshold_10per'] != 0).sum(), size= df_10per.shape[0]))
        return(df_10per)


# Define a function to get list of phenos that has significant kmers associated with
def get_phenos_with_sig_kmers(res_tab, threshold=5):
    res_tab_filt = res_tab.drop(['kmers_pass_threshold_{thr}per'.format(thr=str(threshold))], axis=1)
    res_tab_filt = res_tab_filt[res_tab_filt['num_kmers_pass_threshold_{thr}per'.format(thr=str(threshold))] > 0]
    res_tab_filt.to_csv(args.out_path + "/phenos_with_kmers_pass_threshold_{thr}per.txt".format(thr=str(threshold)), index=False, sep="\t")

    phenos_with_sig_kmers=res_tab_filt["pheno"]
    phenos_with_sig_kmers.to_csv(args.out_path + "/phenos_with_sig_kmers_list_{thr}per.txt".format(thr=str(threshold)), index=False, sep="\t")

## Parse pass_threshold_5per results and write out the output
results_5per = parse_kmersGWAS_results(args.dir_path, threshold=5)
results_5per.to_csv(args.out_path + "/kmers_gwas_results_table_5per.txt", index=False, sep="\t")
## Get phenos with significant k-mers list and write out the output
get_phenos_with_sig_kmers(results_5per, threshold=5)

## If 10 per threshold activated use the results from pass_threshold_10per
if args.use_10_per:
    results_10per = parse_kmersGWAS_results(args.dir_path, threshold=10)
    results_10per.to_csv(args.out_path + "/kmers_gwas_results_table_10per.txt", index=False, sep="\t")
    get_phenos_with_sig_kmers(results_10per, threshold=10)