#!/usr/bin/env python3
from setuptools import setup, find_packages

# Read README file
try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "ONT-optimized circular DNA detection via multi-modal analysis"

setup(
    name="circDNA-detection",  # Note: using hyphen for PyPI name
    version="1.0.0",
    author="Samarth Mathur, PhD",
    author_email="samarth8392@gmail.com",
    description="ONT-optimized circular DNA detection via multi-modal analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samarth8392/circDNA_detection",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pysam>=0.19.0",
        "numpy>=1.19.0",
        "scipy>=1.6.0",
    ],
    entry_points={
        "console_scripts": [
            "circDNA-detect=circDNA_detection.circular_dna_detector:main",
        ],
    },
    include_package_data=True,
)