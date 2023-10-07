import os
from setuptools import setup, find_packages

# Get the long description from the README file
setup_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(setup_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="kgwasflow",
    version="1.3",
    python_requires=">3.10",
    description="kGWASflow is a Snakemake workflow for performing k-mers-based GWAS.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/akcorut/kGWASflow",
    author="Adnan Kivanc Corut",
    keywords="k-mers GWAS genomics snakemake",
    license="MIT",
    packages=find_packages(),
    package_data={'workflow': [
        "Snakefile",
        "scripts/*",
        "envs/*",
        "rules/*",
        "scripts/*",
        "schemas/*",
        "report/*",
        "config/*", 
        "test/config_ecoli/*",
        "test/config_test/*",
        "test/data/ecoli_ref/*", 
        "test/data/ecoli_phenos/*",
        "test/data/test_ref/*",
        "test/data/test_phenos/*",
        "test/data/test_reads/*"]},
    # include_package_data= True,
    entry_points={
        "console_scripts": [
            "kgwasflow = workflow.kgwasflow:main",
        ],
    },
    install_requires=[
        "snakemake==7.25.0",
        "click"
    ],
)