# CircONTrack Peaks - Statistical Coverage Peak Detection

## Overview

The `circontrack-peaks` module provides rigorous statistical peak calling using negative binomial distribution modeling. Unlike simple fold-change thresholds, this module identifies **statistically significant** coverage enrichments with proper multiple testing correction.

## Key Features

- **Negative binomial modeling** - Appropriate for overdispersed sequencing data
- **FDR control** - Benjamini-Hochberg correction for multiple testing
- **Genome-wide analysis** - Process all chromosomes in parallel
- **Peak merging** - Combine nearby significant windows
- **Visualization** - Coverage plots with p-value tracks

## Statistical Framework

### Why Negative Binomial?

ONT sequencing coverage is typically **overdispersed** (variance > mean):

```
Poisson:           Variance = Mean
Negative Binomial: Variance > Mean (overdispersion)

ONT Coverage:      High variance due to:
                   - Read length variation
                   - Mapping biases
                   - Circular DNA amplification
```

### Statistical Test

For each window, we test:

```
H₀: Coverage comes from background distribution
H₁: Coverage is significantly enriched

P-value = P(X ≥ observed | H₀)
```

### Multiple Testing Correction

```python
# Benjamini-Hochberg FDR correction
Adjusted p-value = p-value × (n_tests / rank)

# Keeps false discovery rate < threshold (e.g., 5%)
```

## Usage

### Basic Peak Calling

```bash
# Simple peak detection
circontrack-peaks sample.bam -o peaks.bed

# With reference for accurate chromosome lengths
circontrack-peaks sample.bam -r reference.fa -o peaks.bed
```

### Advanced Parameters

```bash
circontrack-peaks sample.bam -r reference.fa \
  --window-size 500 \          # Smaller windows for higher resolution
  --step-size 250 \            # Overlapping windows
  --fdr 0.01 \                 # Stricter FDR (1%)
  --min-fold 3.0 \             # Require 3-fold enrichment
  --min-coverage 10 \          # Minimum 10X coverage
  --merge-distance 2000 \      # Merge peaks within 2kb
  -o significant_peaks.bed
```

### Specific Chromosomes with Visualization

```bash
# Analyze specific chromosomes and generate plots
circontrack-peaks sample.bam -r reference.fa \
  --chromosomes chr1 chr2 chr3 NC_001234 \
  --plot \
  --plot-dir coverage_plots \
  -o selected_peaks.bed
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--window-size` | 1000 | Window size in bp |
| `--step-size` | 500 | Step size (default: window/2) |
| `--fdr` | 0.05 | False discovery rate threshold |
| `--min-fold` | 2.0 | Minimum fold change over background |
| `--min-coverage` | 5 | Minimum coverage in peak |
| `--merge-distance` | 1000 | Merge peaks within this distance |
| `--plot` | False | Generate coverage plots |

