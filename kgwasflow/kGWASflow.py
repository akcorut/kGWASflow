#!/usr/bin/env python3

import argparse
import os

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run kGWASflow')
    parser.add_argument('-c', '--config-file', default='config/config.yaml', help='Path to the config.yaml file (default: config.yaml)')
    # parser.add_argument('--output', type=str, help='Path to output directory')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--conda-frontend', type=str, default='conda', help='Conda frontend to use')
    parser.add_argument('-n', '--dryrun', action='store_true', help='Dry run', required=False)
    parser.add_argument('--samples', type=str, help='Path to samples.tsv file', required=False)
    parser.add_argument('--phenotypes', type=str, help='Path to phenos.tsv file', required=False)
    parser.add_argument('-r', '--generate-report', action='store_true', help='create an HTML report', required=False)
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity', required=False)
    args = parser.parse_args()

    # if config file is not specified, use the default config file
    if args.config_file is None:
        args.config_file = 'config/config.yaml'

    # Define the command to run snakemake
    cmd = f'snakemake --use-conda --conda-frontend {args.conda_frontend} --cores {args.threads} --rerun-triggers mtime --rerun-incomplete --configfile {args.config_file}'

    if args.dryrun:
        cmd += ' --dryrun'

    # Add the --report flag if specified
    if args.generate_report:
        if args.dryrun:
            cmd += ' --report kGWASflow-report.html'
        if not args.dryrun:
            cmd += ' --dryrun --report kGWASflow-report.html'

    if args.verbose:
        cmd += ' --verbose'
        
    os.system(cmd)
    
if __name__ == "__main__":
    main()