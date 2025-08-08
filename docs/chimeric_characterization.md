# CircONTrack + Characterization Integration Guide

## Understanding the Integration

### File Structure
First, you need to place the characterization module in the CircONTrack directory:

```
CircONTrack/
├── circDNA_detection/
│   ├── __init__.py
│   ├── coverage_analyzer.py
│   ├── junction_detector.py
│   ├── split_read_analyzer.py
│   ├── utils.py
│   ├── confidence_scorer.py
│   └── characterization.py  <-- ADD THIS FILE HERE
└── circontrack.py  <-- Main CircONTrack script
```

## Command-Line Usage Options

### Standard CircONTrack and then Standalone Characterization

Run normal CircONTrack workflow that detects circular DNA:
**Output**: Standard BED file with circular DNA regions detected by coverage/junction/split-read methods

```bash
# Standard CircONTrack command
circontrack input.bam reference.fasta -o circdna_output.bed

# Then characterize the detected regions
python circDNA_detection/characterization.py characterize \
    circdna_regions.bed \
    mm10.fasta \
    mm10_plus_viral.fasta \
    sample_aligned_to_mm10.bam \
    sample_aligned_to_combined.bam \
    --viral-prefix "viral_" \
    -o characterized_results
```
**Output**: 
- `characterized_results_characterized.bed` - Enhanced BED with viral/host classification
- `characterized_results_summary.tsv` - Detailed report

## Complete Workflow Example

### Step 1: Prepare Your Data

```bash
# You need these files:
# 1. Your sequencing data
sample.fastq.gz

# 2. Reference genomes
mm10.fasta                    # Mouse reference only
mm10_plus_viral.fasta         # Combined mouse + viral genomes

# 3. Create alignments to BOTH references
minimap2 -ax map-ont mm10.fasta sample.fastq.gz | samtools sort -o sample.mm10.bam
minimap2 -ax map-ont mm10_plus_viral.fasta sample.fastq.gz | samtools sort -o sample.combined.bam
samtools index sample.mm10.bam
samtools index sample.combined.bam
```

### Step 2: Run Integrated CircONTrack + Characterization

```bash
# With the modified circontrack.py
circontrack sample.combined.bam mm10.fasta \
    --characterize \
    --host-ref mm10.fasta \
    --combined-ref mm10_plus_viral.fasta \
    --host-bam sample.mm10.bam \
    --combined-bam sample.combined.bam \
    --viral-prefix "viral_" \
    --min-fold-enrichment 1.5 \
    --min-coverage 5 \
    -o results
```

## Output Interpretation

The final output tells you:

| Region | Type | Meaning |
|--------|------|---------|
| chr1:1000-5000 | host_only | Normal mouse eccDNA |
| chr2:2000-3000 | viral_only | Pure viral episome |
| chr3:4000-6000 | integration_site | Virus integrated into mouse genome |
| chr4:7000-9000 | chimeric | Complex host-viral rearrangement |

