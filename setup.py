#!/usr/bin/env python3
"""
Setup script for CircONTrack
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="CircONTrack",
    version="0.3.0",  # Bumped version for new features
    author="Samarth Mathur, PhD",
    author_email="samarth8392@gmail.com",
    description="ONT-optimized circular DNA detection with viral classification and assembly preparation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samarth8392/CircONTrack",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pysam>=0.19.0",
        "numpy>=1.19.0",
        "scipy>=1.6.0",
    ],
    extras_require={
        "full": [
            "pandas>=1.0.0",  # For better stats output
            "tqdm>=4.60.0",   # For progress bars
        ],
        "dev": [
            "pytest>=6.0.0",
            "black>=21.0",
            "flake8>=3.9.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "circontrack=circDNA_detection.circular_dna_detector:main",
            "circontrack-classify=circDNA_detection.classify:main",
            "circontrack-assemble=circDNA_detection.assemble:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/samarth8392/CircONTrack/issues",
        "Source": "https://github.com/samarth8392/CircONTrack",
        "Documentation": "https://samarth8392.github.io/CircONTrack/",
    },
)