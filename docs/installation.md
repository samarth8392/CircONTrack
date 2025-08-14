# Installation

This guide covers different methods to install circDNA Detection on your system.

## Requirements

### System Requirements

- **Operating System**: Linux, macOS, or Windows (with WSL)
- **Python**: Python 3.7 or higher
- **Memory**: At least 4GB RAM (8GB+ recommended for large datasets)
- **Storage**: Sufficient space for BAM files and output

### Python Dependencies

The following Python packages are required:

- `pysam >= 0.19.0` - For BAM file handling
- `numpy >= 1.19.0` - For numerical computations
- `scipy >= 1.6.0` - For statistical analysis

## Installation Methods

### Method 1: Install from GitHub (Recommended)

Install directly from the GitHub repository:

```bash
pip install git+https://github.com/samarth8392/CircONTrack.git

This method automatically installs all dependencies and provides the latest stable version.

### Method 2: Development Installation

For development or to contribute to the project:

```bash
# Clone the repository
git clone https://github.com/samarth8392/CircONTrack.git
cd CircONTrack

# Install in development mode
pip install -e .
```

This method allows you to modify the code and have changes reflected immediately.


### Method 2: Conda Environment

If you use Conda, you can create a dedicated environment:

```bash
# Create a new conda environment
conda create -n CircONTrack python=3.8

# Activate the environment
conda activate CircONTrack

# Install dependencies
conda install -c bioconda pysam numpy scipy

# Install circDNA Detection
pip install git+https://github.com/samarth8392/CircONTrack.git
```

## Verification

After installation, verify that CircONTrack is correctly installed

### Command Line Tool

```bash
circontrack -h

    
   _____ _           ____  _   _ _______             _    
  / ____(_)         / __ \| \ | |__   __|           | |   
 | |     _ _ __ ___| |  | |  \| |  | |_ __ __ _  ___| | __
 | |    | | '__/ __| |  | | . ` |  | | '__/ _` |/ __| |/ /
 | |____| | | | (__| |__| | |\  |  | | | | (_| | (__|   < 
  \_____|_|_|  \___|\____/|_| \_|  |_|_|  \__,_|\___|_|\_\ 

    

     Circular DNA Tracking and Detection via ONT-based Multi-modal Signal Integration 

      Version:  0.3.0   
    
usage: circontrack [-h] [-o OUTPUT] [-c CHROMOSOME] [-q] [--log-level {DEBUG,INFO,WARNING,ERROR}] [--min-fold-enrichment MIN_FOLD_ENRICHMENT] [--min-coverage MIN_COVERAGE] [--min-length MIN_LENGTH]
                   [--max-length MAX_LENGTH] [--min-confidence MIN_CONFIDENCE]
                   bam_file reference_file

Detect circular DNA from ONT sequencing data using multi-modal analysis

positional arguments:
  bam_file              Input BAM file
  reference_file        Reference FASTA file

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output BED file
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Analyze specific chromosome only
  -q, --quiet           Disable verbose output (quiet mode)
  --log-level {DEBUG,INFO,WARNING,ERROR}
                        Set logging level
  --min-fold-enrichment MIN_FOLD_ENRICHMENT
                        Minimum fold enrichment for coverage-based detection
  --min-coverage MIN_COVERAGE
                        Minimum coverage threshold
  --min-length MIN_LENGTH
                        Minimum circular DNA length
  --max-length MAX_LENGTH
                        Maximum circular DNA length
  --min-confidence MIN_CONFIDENCE
                        Minimum confidence score threshold

```

### Getting Help

If you encounter issues not covered here:

1. Check the [GitHub Issues](https://github.com/samarth8392/circDNA_detection/issues) page
2. Create a new issue with:
   - Your operating system
   - Python version (`python --version`)
   - Complete error message
   - Installation method used

## Updating

To update to the latest version:

```bash
pip install --upgrade git+https://github.com/samarth8392/CircONTrack.git
```

For development installations:

```bash
cd CircONTrack
git pull origin main
pip install -e .
```

## Uninstalling

To remove circDNA Detection:

```bash
pip uninstall CircONTrack
```