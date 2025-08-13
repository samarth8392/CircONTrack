from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="CircONTrack",
    version="0.2.0",
    author="Samarth Mathur, PhD",
    author_email="samarth8392@gmail.com",
    description="Circular DNA Tracking and Detection via ONT-based Multi-modal Signal Integration",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samarth8392/CircONTrack",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.10",
    install_requires=[
        "pysam>=0.19.0",
        "numpy>=1.19.0",
        "scipy>=1.6.0",
        "tqdm>=4.60.0",  # Added tqdm dependency
        "rich>=10.0.0",  # Added rich dependency for console output
    ],
    entry_points={
        "console_scripts": [
            'circontrack=circDNA_detection.circular_dna_detector:main',
            'circontrack-classify=circDNA_detection.classify:main',
        ],
    },
)