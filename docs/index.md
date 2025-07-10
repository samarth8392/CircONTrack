# circDNA Detection

**ONT-optimized circular DNA detection via multi-modal analysis combining coverage patterns, junction detection, and split-read analysis.**

## Overview

circDNA Detection is a specialized tool designed for identifying circular DNA elements in Oxford Nanopore Technologies (ONT) long-read sequencing data. The package employs a sophisticated multi-modal approach that combines multiple detection strategies to achieve high sensitivity and specificity.

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
pip install git+https://github.com/samarth8392/circDNA_detection.git
```

### Basic Usage

```bash
circDNA-detect input.bam reference.fasta -o output.bed
```

### Python API

```python
from circDNA_detection import CircularDNADetector

detector = CircularDNADetector(
    min_fold_enrichment=1.5,
    min_coverage=5,
    min_length=200
)
candidates = detector.detect_circular_dna("input.bam", "reference.fasta")
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