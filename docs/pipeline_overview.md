# CircDNA Detection Package Analysis

## Overview

The `circontrack` command is a Python workflow for detecting candidate circular DNA intervals in Oxford Nanopore Technologies (ONT) long-read sequencing data. The package combines three implemented evidence streams: coverage patterns, junction-like signatures, and split-read signatures. This repository does not provide benchmarked sensitivity or specificity estimates.

## Package Architecture

### Core Components

The package implements a comprehensive detection pipeline that consists of four main phases:

1. **Coverage Pattern Analysis**
2. **Junction Detection** 
3. **Split-Read Analysis**
4. **Multi-Modal Integration**

## Detailed Function Analysis

### 1. Coverage Pattern Analysis

**Purpose**: Identifies regions with elevated coverage that may indicate circular DNA amplification.

**Logic**: 
- Calculates coverage depth across genomic regions
- Identifies regions with coverage significantly higher than background
- Uses configurable fold-enrichment thresholds (default: 1.5x)
- Filters based on minimum coverage depth (default: 5x)

**Uses median and MAD (Median Absolute Deviation) instead of mean/standard deviation for robustness against outliers**

**Implementation note**:

This approach tests for elevated read counts relative to local/background windows. Elevated coverage can be consistent with circular DNA or amplification, but it can also arise from mapping and repetitive-sequence artifacts.

### 2. Junction Detection

**Purpose**: Identifies back-to-back junction signatures characteristic of circular DNA.

**Logic**:
- Searches for reads that span the junction point where the circular DNA "loops back"
- Detects characteristic back-to-back alignments
- Validates junction signatures through read orientation analysis

**Implementation note**:
Junction-like read signatures provide direct evidence for a breakpoint pattern consistent with circular topology. The current implementation uses SA tags, soft clips, and large CIGAR gaps/skips; it does not realign clipped sequence.

### 3. Split-Read Analysis

**Purpose**: Analyzes split alignments to identify circular DNA signatures.

**Logic**:
- Examines reads that align to multiple locations
- Identifies split alignments that suggest circular topology
- Validates split-read patterns consistent with circular DNA structure

**Implementation note**:
Split-read analysis uses SA-tagged supplementary alignments. These reads can support circular-junction hypotheses, but interpretation depends on alignment quality and repeat context.

### 4. Multi-Modal Integration

**Purpose**: Combines evidence from all three detection methods and generates confidence scores.

**Logic**:
- Integrates results from coverage, junction, and split-read analyses
- Assigns confidence scores based on multiple evidence types
- Filters candidates based on configurable thresholds
- Outputs results in standard BED format with additional annotation

**Implementation note**:
The integrator merges nearby candidates and scores available evidence. A multi-method candidate has more implemented evidence fields, but the score is not a calibrated probability.

## Configuration Parameters

### Key Parameters and Their Logic

| Parameter | Default | Purpose | Logic Assessment |
|-----------|---------|---------|------------------|
| `min_fold_enrichment` | 1.5 | Minimum coverage fold increase | Used by coverage detection thresholding |
| `min_coverage` | 5 | Minimum coverage threshold | Used by coverage detection thresholding |
| `min_length` | 200 | Minimum circular DNA length | Accepted by CLI; final orchestration does not currently apply an extra length filter |
| `max_length` | 100,000 | Maximum circular DNA length | Accepted by CLI; module-specific hard-coded checks still apply |


## Output Format

The package outputs results in BED format with additional columns:
- Standard BED columns (chr, start, end, name, score, strand)
- Detection method information
- Confidence scores
- Additional details

## Strengths of the Implementation

1. **ONT-Oriented Evidence**: Uses signals commonly available in long-read alignments
2. **Multi-Modal Approach**: Records coverage, junction, and split-read evidence types
3. **Configurable Thresholds**: Allows adaptation to different experimental conditions
4. **Comprehensive Scoring**: Provides confidence measures for downstream analysis
5. **Standard Output**: Uses widely-accepted BED format for compatibility

## Potential Considerations

1. **Parameter Sensitivity**: Candidate ranking and filtering depend on dataset-specific thresholds
2. **Computational Complexity**: Multi-modal analysis may be computationally intensive for large datasets
3. **Artifact Rate**: Repeats, mapping artifacts, and coverage nonuniformity can produce candidate-like signals

## Dependencies and Requirements

- Python ≥ 3.10
- pysam ≥ 0.19.0 (for BAM/SAM file handling)
- numpy ≥ 1.19.0 (for numerical computations)
- scipy ≥ 1.6.0 (for statistical analysis)
