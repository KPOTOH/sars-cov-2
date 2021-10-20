from setuptools import find_packages, setup

setup(
    name="sars_processing",
    packages=find_packages(),
    version="0.1.0",
    description="Process sars-cov-2 sequences from gisaid and get mutations from builded phylogenetic tree",
    author="Bogdan Efimenko",
    install_requires=[
        "matplotlib==3.4.3",
        "numpy==1.21.2",
        "pandas==1.3.4",
        "biopython==1.79",
        "tqdm==4.62.3",
        "ete3==3.1.2",
        "click==8.0.3",
    ],
    license="MIT",
)