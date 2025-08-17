<p align="center">
  <img src="./docs/img/CircONTrack_logo.png" width="1000">
</p>

# CircONTrack

**ONT-optimized circular DNA detection via multi-modal analysis combining coverage patterns, junction detection, and split-read analysis.**

## Overview

CircONTrack is a specialized tool designed for identifying circular DNA elements in Oxford Nanopore Technologies (ONT) long-read sequencing data. The package employs a sophisticated multi-modal approach that combines multiple detection strategies to achieve high sensitivity and specificity.

## Key Features

- **Multi-modal detection**: Combines coverage analysis, junction detection, and split-read analysis
- **ONT-optimized**: Specifically designed for Oxford Nanopore Technologies long-read sequencing
- **High sensitivity**: Detects circular DNA elements with configurable thresholds  
- **Comprehensive scoring**: Multi-evidence confidence scoring system
- **Flexible output**: Results in standard BED format with additional metadata

## Detection Methods

The pipeline integrates four complementary detection approaches:

1. **Coverage Pattern Analysis** - Identifies regions with elevated coverage patterns characteristic of circular DNA
2. **Junction Detection** - Finds back-to-back junction signatures at breakpoints
3. **Split-Read Analysis** - Analyzes split alignments for circular signatures
4. **Multi-Modal Integration** - Combines evidence from all methods and scores candidates

## Quick Start

### Installation

```bash
pip install git+https://github.com/samarth8392/CircONTrack.git
```

### Basic Usage

```bash
circontrack input.bam reference.fasta -o output.bed

circontrack -h

    
   _____ _           ____  _   _ _______             _    
  / ____(_)         / __ \| \ | |__   __|           | |   
 | |     _ _ __ ___| |  | |  \| |  | |_ __ __ _  ___| | __
 | |    | | '__/ __| |  | | . ` |  | | '__/ _` |/ __| |/ /
 | |____| | | | (__| |__| | |\  |  | | | | (_| | (__|   < 
  \_____|_|_|  \___|\____/|_| \_|  |_|_|  \__,_|\___|_|\_\ 

    

     Circular DNA Tracking and Detection via ONT-based Multi-modal Signal Integration 

      Version: 0.4.0 
    
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

## Requirements

- Python ≥ 3.7
- pysam ≥ 0.19.0
- numpy ≥ 1.19.0
- scipy ≥ 1.6.0

## Support

For questions, issues, or contributions, please visit our [GitHub repository](https://github.com/samarth8392/circDNA_detection) or submit an issue on the [issue tracker](https://github.com/samarth8392/circDNA_detection/issues).

## Citation

If you use circDNA Detection in your research, please cite:

```
[Citation information will be added upon publication]
```

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/samarth8392/circDNA_detection/blob/main/LICENSE) file for details.