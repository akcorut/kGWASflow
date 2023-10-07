import sys
import os
import subprocess
import click
import logging

# Get the directory path of this file
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
workflow_dir = os.path.join(base_dir, "workflow")

logging.basicConfig(level=logging.INFO)

def get_snakefile(file="Snakefile"):
    snake_file = os.path.join(workflow_dir, file)
    if not os.path.exists(snake_file):
        sys.exit("Unable to locate the  Snakefile;  tried %s" % snake_file)
    return snake_file

def get_configfile(file="config.yaml"):
    config_file = os.path.join(workflow_dir, "config", file)
    if not os.path.exists(config_file):
        sys.exit("Unable to locate the config.yaml file;  tried %s" % config_file)
    return config_file

def get_samplefile(file="samples.tsv"):
    sample_file = os.path.join(workflow_dir, "config", file)
    if not os.path.exists(sample_file):
        sys.exit("Unable to locate the samples.tsv file;  tried %s" % sample_file)
    return sample_file

def get_phenosfile(file="phenos.tsv"):
    phenos_file = os.path.join(workflow_dir, "config", file)
    if not os.path.exists(phenos_file):
        sys.exit("Unable to locate the phenos.tsv file;  tried %s" % phenos_file)
    return phenos_file

def get_testdir(dir="test"):
    test_dir = os.path.join(workflow_dir, dir)
    if not os.path.exists(test_dir):
        sys.exit("Unable to locate the test directory;  tried %s" % test_dir)
    return test_dir

def show_help_message():
    message = (
        "\nCLUSTER EXECUTION:\n"
        "\n kgwasflow run ... --profile [profile]\n"
        "\n For information on Snakemake profiles see:\n"
        "https://snakemake.readthedocs.io/en/stable/executing/cluster.html"
        "\nhttps://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles\n"
        "\nUSAGE EXAMPLES:\n"
        "\n    kgwasflow init [OPTIONS]        Initialize a new kGWASflow working directory\n"
        "\n    kgwasflow run [OPTIONS]         Run the kGWASflow workflow\n"
        "\n    kgwasflow test [OPTIONS]        Run the kGWASflow test\n"
        "\n    kgwasflow --help"
        "\n\nInit examples:"
        "\n\n1. Initialize a new kGWASflow working directory in the current directory:"
        "\n\n    kgwasflow init"
        "\n\n2. Initialize a new kGWASflow working directory in a specified directory:"
        "\n\n    kgwasflow init --work-dir path/to/working_dir"
        "\n\nRun examples:"
        "\n\n1. Run kGWASflow with the default config file, default arguments and 16 threads:\n"
        "\n    kgwasflow run -t 16 --snake-default"
        "\n\n2. Run kGWASflow with a custom config file and default settings:\n"
        "\n    kgwasflow run -t 16 -c path/to/custom_config.yaml"
        "\n\n3. Run kGWASflow with user defined output directory:\n"
        "\n    kgwasflow run -t 16 --output path/to/output_dir"
        "\n\n4. Run kGWASflow in dryrun mode to see what tasks would be executed:\n"
        "\n    kgwasflow run -t 16 -n"
        "\n\n5. Run kGWASflow using mamba as the conda frontend:\n"
        "\n    kgwasflow run -t 16 --conda-frontend mamba"
        "\n\n6. Run kGWASflow and generate an HTML report:\n"
        "\n    kgwasflow run -t 16 --generate-report"
        "\n\nTest examples:"
        "\n\n1. Run the kGWASflow test in dryrun mode to see what tasks would be executed:\n"
        "\n    kgwasflow test -t 16 -n"
        "\n\n2. Run the kGWASflow test using the default 'test' dataset with 16 threads:\n"
        "\n    kgwasflow test -t 16"
        "\n\n3. Run the kGWASflow test and define the working directory:\n"
        "\n    kgwasflow test -t 16 --work-dir path/to/test_output_dir"
        "\n\n4. Run the kGWASflow test using the 'ecoli' dataset:\n"
        "\n    kgwasflow test -t 16 --dataset ecoli"
        "\n\n5. Run the kGWASflow test using the default 'test' dataset:\n"
        "\n    kgwasflow test -t 16 --dataset test"
    )
    return message



def show_ascii_art():
    click.echo("""
    \b           
     _     _______          __      _____  __ _               
    | |   / ____\ \        / /\    / ____|/ _| |              
    | | _| |  __ \ \  /\  / /  \  | (___ | |_| | _____      __
    | |/ / | |_ | \ \/  \/ / /\ \  \___ \|  _| |/ _ \ \ /\ / /
    |   <| |__| |  \  /\  / ____ \ ____) | | | | (_) \ V  V / 
    |_|\_\\_____|   \/  \/_/    \_\_____/|_| |_|\___/ \_/\_/  
    \b
    kGWASflow: A Snakemake Workflow for k-mers Based GWAS
    """)

def run_snake(snakefile, config_file, threads, work_dir, conda_frontend, dryrun, generate_report, snake_default, rerun_triggers, verbose, unlock, snakemake_args):
    # Define the command to run snakemake
    cmd = ['snakemake', '--use-conda', '--conda-frontend', conda_frontend, '--cores', str(threads)]

    if snakefile:
        cmd += ['--snakefile', snakefile]
    
    # if config file is provided, use it
    if config_file:
        cmd += ['--configfile', config_file]
        
    # if output directory is provided, use it
    if work_dir:
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        cmd += ['--directory', work_dir]
        
    if dryrun:
        cmd.append('--dryrun')

    if generate_report:
        if dryrun:
            cmd.append('--report')
            cmd.append('kGWASflow-report.html')
        if not dryrun:
            cmd.append('--dryrun')
            cmd.append('--report')
            cmd.append('kGWASflow-report.html')
    
    if snakemake_args:
        cmd += snakemake_args
    
    if rerun_triggers:
        cmd += ['--rerun-triggers'] + list(rerun_triggers)
    
    if snake_default:
        default_snakemake_args = ["--rerun-incomplete", "--printshellcmds", "--nolock", "--show-failed-logs"]
        cmd += default_snakemake_args
    
    if verbose:
        cmd.append('--verbose')
        
    if unlock:
        cmd.append('--unlock')
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("Error running Snakemake: {}".format(e))