import os
from glob import glob
from setuptools import setup

exec(open("synspace/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="synspace",
    version=__version__,
    description="Generative model based on forward synthesis rules",
    author="Andrew White",
    author_email="andrew.white@rochester.edu",
    url="https://github.com/whitead/synspace",
    license="MIT",
    packages=["synspace", "synspace/reos"],
    install_requires=[
        "rdkit",
        "tqdm",
        "numpy",
        "importlib_resources",
        "requests",
        "pandas",
        "pystow",
        "skunk",
    ],
    package_data={"synspace": ["rxn_data/*.json", "rxn_data/*.bz2"]},
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
