# circDNA Detection

ONT-optimized circular DNA detection via multi-modal analysis combining coverage patterns, junction detection, and split-read analysis.

## Features

- **Multi-modal detection**: Combines coverage analysis, junction detection, and split-read analysis
- **ONT-optimized**: Specifically designed for Oxford Nanopore Technologies long-read sequencing
- **High sensitivity**: Detects circular DNA elements with configurable thresholds
- **Comprehensive scoring**: Multi-evidence confidence scoring system

## Installation

### From GitHub (Recommended)

```bash
pip install git+https://github.com/samarth8392/circDNA_detection.git
```

### From Source

```bash
git clone https://github.com/samarth8392/circDNA_detection.git
cd circDNA_detection
pip install -e .
```

## Requirements

- Python ≥ 3.7
- pysam ≥ 0.19.0
- numpy ≥ 1.19.0
- scipy ≥ 1.6.0

## Usage

### Command Line

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

## Parameters

- `--min-fold-enrichment`: Minimum fold enrichment for coverage detection (default: 1.5)
- `--min-coverage`: Minimum coverage depth (default: 5)
- `--min-length`: Minimum circular DNA length (default: 200)
- `--max-length`: Maximum circular DNA length (default: 100000)
- `-c, --chromosome`: Analyze specific chromosome only
- `-o, --output`: Output BED file (default: circular_dna_ont.bed)

## Output Format

Results are written in BED format with additional columns:
- Standard BED columns (chr, start, end, name, score, strand)
- Detection method information
- Confidence scores
- Additional details

## Algorithm

The detection pipeline consists of four phases:

1. **Coverage Pattern Analysis**: Identifies regions with elevated coverage
2. **Junction Detection**: Finds back-to-back junction signatures
3. **Split-Read Analysis**: Analyzes split alignments for circular signatures
4. **Multi-Modal Integration**: Combines evidence and scores candidates


## Issues

Please report issues on the [GitHub issue tracker](https://github.com/samarth8392/circDNA_detection/issues).

