#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=kGWASflow
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=500gb

#SBATCH --mail-user=ac32082@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=kGWASflow.%j.out
#SBATCH --error=kGWASflow.%j.err
#SBATCH --export=NONE

source /apps/lmod/lmod/init/zsh

cd /scratch/ac32082/kGWASflow/workflow

module load Anaconda3/5.0.1
source activate kGWASflow

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --report report.html --reason --use-conda --conda-frontend conda --verbose --cores 30 --rerun-incomplete -s Snakefile