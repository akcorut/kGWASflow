#!/usr/bin/env python3

import click
import os
import shutil
from shutil import copyfile

from workflow import __version__
from .cli_utils import run_snake, get_phenosfile, get_samplefile, get_testdir ,get_snakefile, get_configfile, show_ascii_art, show_help_message, workflow_dir, base_dir

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
def cli():
    """kGWASflow is a Snakemake workflow for performing k-mers-based GWAS.

    For more options, run:
    kgwasflow --help
    kgwasflow init --help
    kgwasflow run --help
    kgwasflow test --help"""
    pass

def common_options(func):
    """Decorator for common options."""
    options = [
        click.option('-s', '--snakefile', type=click.Path(dir_okay=True, writable=True, resolve_path=True), help='Path to the Snakefile.'),
        click.option('-c', '--config-file', type=click.Path(dir_okay=True, writable=True, resolve_path=True), help='Path to the config.yaml file'),
        click.option('-t', '--threads', default=8, type=int, help='Number of threads (default: 8).', show_default=True),
        click.option('-d', '--work-dir', help="kGWASflow working directory.", type=click.Path(dir_okay=True, writable=True, readable=True)),
        click.option('--conda-frontend', default='conda', type=str, help='Conda frontend to use.'),
        click.option('-n', '--dryrun', is_flag=True, default=False, show_default=True, help='Dry run. Do not execute the workflow, but show which jobs would be executed.'),
        click.option('-r', '--generate-report', is_flag=True, default=False, help='Create a kGWASflow HTML report.', show_default=True),
        click.option('--snake-default', is_flag=True, default=False, help='Useful default snakemake arguments.', show_default=True),
        click.option('--rerun-triggers', multiple=True, default= ["mtime", "params", "input", "software-env", "code"], help='Rerun all jobs that have at least one of the specified trigger files changed.', show_default=True),
        click.option('--unlock', is_flag=True, help='Unlock the workflow if it is locked.'),
        click.option('-v', '--verbose', is_flag=True, default=False, help='Increase output verbosity.'),
        click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED),
    ]
    for option in reversed(options):
        func = option(func)
    return func

@click.command(epilog=show_help_message(), context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@common_options
def run(snakefile, config_file, **kwargs):
    """Run kGWASflow workflow."""
    if not snakefile:
        snakefile = get_snakefile()
    run_snake(snakefile, config_file, **kwargs)
    
@common_options
@click.command(epilog=show_help_message(), context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
def test(snakefile, config_file, **kwargs):
    """Test kGWASflow workflow."""
    if not snakefile:
        snakefile = get_snakefile()
    
    test_config_file = os.path.join(workflow_dir, "test", "config_ecoli", "config.yaml")
    if not config_file:
        config_file = test_config_file
    run_snake(snakefile, config_file, **kwargs)
    
@cli.command('init', 
             help="Initialize a new kGWASflow working directory with default configuration files.", 
             short_help="Initialize a new kGWASflow working directory."
)
@click.option('-d', '--work-dir', 
              type=click.Path(dir_okay=True, writable=True, resolve_path=True),
              help="Path to the new working directory. If not specified, the current directory will be used.",
              default="."
)
def init_workdir(work_dir):
    """Initialize a new kGWASflow working directory with default configuration files."""
    
    config_file = get_configfile() # Get the default config.yaml file
    sample_file = get_samplefile() # Get the default samples.tsv file
    phenos_file = get_phenosfile() # Get the default phenos.tsv file
    test_dir = get_testdir() # Get the test directory
    
    # Use f-strings to construct paths
    config_path = f"{work_dir}/config"
    config_yaml_path = f"{config_path}/config.yaml"
    samples_tsv_path = f"{config_path}/samples.tsv"
    phenos_tsv_path = f"{config_path}/phenos.tsv"
    new_test_dir = f"{work_dir}/{os.path.basename(test_dir)}"
    
    # Create the new working directory and copy the default files
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(config_path):
        os.makedirs(config_path)
    
    copyfile(config_file, config_yaml_path)
    copyfile(sample_file, samples_tsv_path)
    copyfile(phenos_file, phenos_tsv_path)
    shutil.copytree(test_dir, new_test_dir, dirs_exist_ok=True)
    
cli.add_command(run)
cli.add_command(test)

def main():
    show_ascii_art()
    cli()
    
if __name__ == "__main__":
    main()