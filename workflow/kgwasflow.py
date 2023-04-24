#!/usr/bin/env python3

import click
import os

from workflow import __version__
from .cli_utils import run_snake, get_snakefile, get_configfile, show_ascii_art, show_help_message, workflow_dir, base_dir

# kgwasflow_dir = os.path.join(os.getcwd())

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
def cli():
    """kGWASflow is a Snakemake workflow for performing k-mers-based GWAS.
    \b
    For more options, run:
    kgwasflow --help
    kgwasflow run --help
    or,
    kgwasflow test --help"""
    pass

def common_options(func):
    """Decorator for common options."""
    options = [
        click.option('-s', '--snakefile', type=click.Path(dir_okay=True, writable=True, resolve_path=True), help='Path to the Snakefile.'),
        click.option('-c', '--config-file', type=click.Path(dir_okay=True, writable=True, resolve_path=True), help='Path to the config.yaml file'),
        click.option('-t', '--threads', default=8, type=int, help='Number of threads (default: 8).', show_default=True),
        click.option('-o', '--output', help="Output directory.", type=click.Path(dir_okay=True, writable=True, readable=True)),
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
    if not config_file:
        config_file = get_configfile()
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
    
cli.add_command(run)
cli.add_command(test)

def main():
    show_ascii_art()
    cli()
    
if __name__ == "__main__":
    main()