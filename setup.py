import os
from glob import glob
from setuptools import setup

exec(open("syngen/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="syngen",
    version=__version__,
    description="Generative model based on forward synthesis rules",
    author="Andrew White",
    author_email="andrew.white@rochester.edu",
    url="https://github.com/whitead/syngen",
    license="MIT",
    packages=["syngen"],
    install_requires=["rdkit", "tqdm", "molbloom",
    "numpy", "importlib_resources", "requests", "useful_rdkit_utils"],
    package_data={"syngen": ["rxn_data/*.json", "rxn_data/*.bz2"]},
    test_suite="tests",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
1