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
pip install git+https://github.com/samarth8392/circDNA_detection.git
```

This method automatically installs all dependencies and provides the latest stable version.

### Method 2: Development Installation

For development or to contribute to the project:

```bash
# Clone the repository
git clone https://github.com/samarth8392/circDNA_detection.git
cd circDNA_detection

# Install in development mode
pip install -e .
```

This method allows you to modify the code and have changes reflected immediately.

### Method 3: Virtual Environment Installation

It's recommended to install in a virtual environment to avoid dependency conflicts:

```bash
# Create a virtual environment
python -m venv circDNA_env

# Activate the virtual environment
# On Linux/macOS:
source circDNA_env/bin/activate
# On Windows:
circDNA_env\Scripts\activate

# Install circDNA Detection
pip install git+https://github.com/samarth8392/circDNA_detection.git
```

### Method 4: Conda Environment

If you use Conda, you can create a dedicated environment:

```bash
# Create a new conda environment
conda create -n circDNA python=3.8

# Activate the environment
conda activate circDNA

# Install dependencies
conda install -c bioconda pysam numpy scipy

# Install circDNA Detection
pip install git+https://github.com/samarth8392/circDNA_detection.git
```

## Verification

After installation, verify that circDNA Detection is correctly installed:

### Command Line Tool

```bash
circDNA-detect --help
```

You should see the help message with available options.

### Python Import

```python
import circDNA_detection
print(circDNA_detection.__version__)
```

### Quick Test

Create a simple test to ensure everything works:

```python
from circDNA_detection import CircularDNADetector

# Create detector instance
detector = CircularDNADetector()
print("Installation successful!")
```

## Troubleshooting

### Common Issues

#### Issue: `ModuleNotFoundError: No module named 'pysam'`

**Solution**: Install pysam manually:
```bash
pip install pysam
```

#### Issue: `ImportError: libhts.so.3: cannot open shared object file`

**Solution**: Install system dependencies:
```bash
# On Ubuntu/Debian:
sudo apt-get install libhts-dev

# On CentOS/RHEL:
sudo yum install htslib-devel

# On macOS:
brew install htslib
```

#### Issue: Permission denied during installation

**Solution**: Use `--user` flag:
```bash
pip install --user git+https://github.com/samarth8392/circDNA_detection.git
```

#### Issue: Version conflicts

**Solution**: Use a fresh virtual environment:
```bash
python -m venv fresh_env
source fresh_env/bin/activate  # On Linux/macOS
pip install git+https://github.com/samarth8392/circDNA_detection.git
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
pip install --upgrade git+https://github.com/samarth8392/circDNA_detection.git
```

For development installations:

```bash
cd circDNA_detection
git pull origin main
pip install -e .
```

## Uninstalling

To remove circDNA Detection:

```bash
pip uninstall circDNA-detection
```