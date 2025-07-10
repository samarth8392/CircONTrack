# Detection Workflow

This page describes the complete workflow of the circDNA detection pipeline, showing how the different detection methods work together to identify circular DNA elements.

## Pipeline Overview

The circDNA detection pipeline consists of four main phases that work together to provide comprehensive circular DNA detection:

```mermaid
graph TD
    A[Input BAM File] --> B[Coverage Pattern Analysis]
    A --> C[Junction Detection]
    A --> D[Split-Read Analysis]
    
    B --> E[Coverage Candidates]
    C --> F[Junction Candidates]
    D --> G[Split-Read Candidates]
    
    E --> H[Multi-Modal Integration]
    F --> H
    G --> H
    
    H --> I[Confidence Scoring]
    I --> J[Output Filtering]
    J --> K[BED Format Output]
    
    style A fill:#e1f5fe
    style K fill:#e8f5e8
    style H fill:#fff3e0
```

## Detailed Workflow

### Phase 1: Coverage Pattern Analysis

```mermaid
graph LR
    A[BAM Input] --> B[Calculate Coverage Depth]
    B --> C[Identify Enriched Regions]
    C --> D[Apply Fold Enrichment Threshold]
    D --> E[Coverage Candidates]
    
    style A fill:#e1f5fe
    style E fill:#e8f5e8
```

**Purpose**: Identifies regions with elevated coverage patterns characteristic of circular DNA elements.

**Process**:
1. Calculate coverage depth across all genomic regions
2. Identify regions with significantly elevated coverage
3. Apply minimum fold enrichment threshold (default: 1.5x)
4. Filter by minimum coverage depth (default: 5x)

### Phase 2: Junction Detection

```mermaid
graph LR
    A[BAM Input] --> B[Identify Split Alignments]
    B --> C[Find Back-to-Back Junctions]
    C --> D[Validate Junction Signatures]
    D --> E[Junction Candidates]
    
    style A fill:#e1f5fe
    style E fill:#e8f5e8
```

**Purpose**: Detects back-to-back junction signatures at circular DNA breakpoints.

**Process**:
1. Identify reads with split alignments
2. Look for back-to-back junction patterns
3. Validate junction signatures for circular characteristics
4. Filter by junction quality and support

### Phase 3: Split-Read Analysis

```mermaid
graph LR
    A[BAM Input] --> B[Extract Split Reads]
    B --> C[Analyze Alignment Patterns]
    C --> D[Identify Circular Signatures]
    D --> E[Split-Read Candidates]
    
    style A fill:#e1f5fe
    style E fill:#e8f5e8
```

**Purpose**: Analyzes split alignments for signatures consistent with circular DNA structures.

**Process**:
1. Extract reads with split alignments
2. Analyze alignment patterns for circular signatures
3. Identify reads supporting circular structures
4. Validate split-read evidence

### Phase 4: Multi-Modal Integration

```mermaid
graph TD
    A[Coverage Candidates] --> D[Evidence Integration]
    B[Junction Candidates] --> D
    C[Split-Read Candidates] --> D
    
    D --> E[Confidence Scoring]
    E --> F[Threshold Filtering]
    F --> G[Final Candidates]
    
    style D fill:#fff3e0
    style G fill:#e8f5e8
```

**Purpose**: Combines evidence from all detection methods and assigns confidence scores.

**Process**:
1. Integrate candidates from all detection methods
2. Calculate multi-evidence confidence scores
3. Apply final filtering thresholds
4. Generate final candidate list

## Parameter Configuration

The workflow can be customized through various parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_fold_enrichment` | 1.5 | Minimum fold enrichment for coverage detection |
| `min_coverage` | 5 | Minimum coverage depth |
| `min_length` | 200 | Minimum circular DNA length |
| `max_length` | 100000 | Maximum circular DNA length |

## Output Generation

The final step generates results in BED format with additional columns:

- Standard BED columns (chr, start, end, name, score, strand)
- Detection method information
- Confidence scores
- Supporting evidence details

## Quality Control

Throughout the pipeline, multiple quality control measures ensure reliable detection:

- **Coverage validation**: Ensures sufficient read support
- **Junction validation**: Confirms genuine circular junction signatures
- **Split-read validation**: Validates split-read evidence quality
- **Multi-modal consensus**: Requires evidence from multiple detection methods for high-confidence calls