import csv
import os, argparse, glob, shutil
import pandas as pd
import subprocess 

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--in_path', help='Input path for reads with k-mers', type=str, required=True)
parser.add_argument('-p', '--pheno', help='Pnenotype name', type=str, required=True)
parser.add_argument('-r', '--ref_index', help='Index path for the reference genome', type=str, required=True)
parser.add_argument('-o', '--out_dir', help='Output directory path', type=str, required=True)
parser.add_argument('-t', '--threads', help='Number of threads to use', type=int, required=True)
args = parser.parse_args()

        ##############################################################
        ############ Align reads with significant k-mers #############
        ##############################################################

# Fetching paths for reads (R1 and R2) with k-mers for a given phenotype
R1_list = glob.glob(args.in_path + '/' + args.pheno + '/' + "*_R1.fastq")

R2_list = glob.glob(args.in_path + '/' + args.pheno + '/' + "*_R2.fastq")

# Generating list of samples that have reads with k-mers for a given phenotype
accession_list=[]
for fastq in R1_list:
    fname = os.path.basename(fastq)
    accession_list.append(fname.split('_')[0])

# Aligning reads for each sample
for acc in accession_list:
    r1_path= args.in_path + '/' + args.pheno + '/' + acc + "_reads_with_kmers_R1.fastq"
    r2_path= args.in_path + '/' + args.pheno + '/' + acc + "_reads_with_kmers_R2.fastq"

    # Making output directory if doesn't exist
    if not os.path.exists(args.out_dir + '/' + args.pheno):
        os.makedirs(args.out_dir + '/' + args.pheno)

    print("################## Aligning reads with k-mers using bowtie2... ##################")
    print("")
    print("Sample: " + acc)
    print("")
    print("R1: " + r1_path)
    print("R2: " + r2_path)
    print("")
    print("")
    subprocess.run(" bowtie2 -p {threads} --very-sensitive-local -x {index} -1 {r1} -2 {r2} -S {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.sam".format(
        r1= r1_path, r2= r2_path, index= args.ref_index, 
        pheno=args.pheno, outdir=args.out_dir, 
        accession= acc, threads=args.threads), shell=True)

    print("################## Converting SAM to BAM... ##################")
    print("")
    print("Sample: " + acc)
    print("")
    subprocess.run(" samtools view -S -b {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.sam > {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.bam".format(
        pheno=args.pheno, outdir=args.out_dir, accession= acc), shell=True)
    print("")

    print("################## Sorting alignment BAM... ##################")
    print("")
    print("Sample: " + acc)
    print("")
    subprocess.run(" samtools sort {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.bam -o {outdir}/{pheno}/{accession}_reads_with_kmers_alignment_sorted.bam".format(
        pheno=args.pheno, outdir=args.out_dir, accession= acc), shell=True)
    print("")    

    print("################## Indexing alignment BAM... ##################")
    print("")
    print("Sample: " + acc)
    print("")
    subprocess.run(" samtools index {outdir}/{pheno}/{accession}_reads_with_kmers_alignment_sorted.bam".format(
        pheno=args.pheno, outdir=args.out_dir, accession= acc), shell=True)
    print("")

    print("################## Removing alignment SAM... ##################")
    print("")
    print("Sample: " + acc)
    print("")   
    subprocess.run(" rm {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.sam".format(
        pheno=args.pheno, outdir=args.out_dir, accession= acc), shell=True)

    subprocess.run(" rm {outdir}/{pheno}/{accession}_reads_with_kmers_alignment.bam".format(
        pheno=args.pheno, outdir=args.out_dir, accession= acc), shell=True)
    print("")
    print("Aligning reads with k-mers is done for the sample: " + acc)

print("")
print("All alignments DONE.")