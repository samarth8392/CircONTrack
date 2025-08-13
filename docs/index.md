# CircONTrack

**ONT-optimized circular DNA detection with viral/host classification via multi-modal analysis combining coverage patterns, junction detection, and split-read analysis.**

## Overview

CircONTrack is a specialized tool for identifying and classifying circular DNA elements in Oxford Nanopore Technologies (ONT) long-read sequencing data. The package employs a sophisticated multi-modal approach to detect circular DNA and can distinguish between host eccDNA, viral episomes, and integration sites.

## Key Features

- **Multi-modal detection**: Combines coverage analysis, junction detection, and split-read analysis
- **ONT-optimized**: Specifically designed for long-read sequencing characteristics
- **Viral/host classification**: Identifies viral episomes, integration sites, and host eccDNA
- **High sensitivity**: Detects circular DNA elements with configurable thresholds  
- **Comprehensive scoring**: Multi-evidence confidence scoring system
- **Standard output**: Results in BED format with detailed annotations

## Installation

```bash
# Install from GitHub
pip install git+https://github.com/samarth8392/CircONTrack.git

# Or clone and install locally
git clone https://github.com/samarth8392/CircONTrack.git
cd CircONTrack
pip install -e .
```

## Quick Start Guide

### Step 1: Prepare Your Reference

For viral detection, combine your host and viral references:

```bash
# Combine host and viral genomes
cat mouse_mm10.fa viral_sequences.fa > mm10_plus_viral.fa
samtools faidx mm10_plus_viral.fa
```

### Step 2: Align Your ONT Reads

```bash
# Align to combined reference
minimap2 -ax map-ont mm10_plus_viral.fa reads.fastq.gz | \
  samtools sort -o sample.bam
samtools index sample.bam
```

### Step 3: Detect Circular DNA

```bash
# Basic detection with default parameters
circontrack sample.bam mm10_plus_viral.fa -o circdna.bed

# Recommended parameters for comprehensive detection
circontrack sample.bam mm10_plus_viral.fa \
  --min-fold-enrichment 2.0 \
  --min-coverage 10 \
  --min-length 200 \
  --max-length 500000 \
  -o circdna.bed
```

### Step 4: Classify Detected Regions (Optional)

If you have viral sequences in your reference, classify the detected regions:

```bash
# Classify as host eccDNA, viral episome, or integration site
circontrack-classify circdna.bed mm10_plus_viral.fa sample.bam \
  -o classified_results
```

## Detailed Usage

### Detection Parameters

```bash
circontrack sample.bam reference.fa \
  --min-fold-enrichment 2.0 \   # Fold enrichment over background (default: 1.5)
  --min-coverage 10 \            # Minimum read coverage (default: 5)
  --min-length 200 \             # Minimum circle size in bp (default: 200)
  --max-length 500000 \          # Maximum circle size in bp (default: 100000)
  --min-confidence 0.3 \         # Minimum confidence score (default: 0.3)
  -c chr19 \                     # Analyze specific chromosome only (optional)
  -o output.bed                  # Output file
```

### Parameter Guidelines

| Target Type | min-fold-enrichment | min-coverage | min-length | Use Case |
|------------|-------------------|--------------|------------|----------|
| **High-confidence eccDNA** | 5.0 | 20 | 1000 | Conservative, low false positives |
| **Comprehensive detection** | 2.0 | 10 | 200 | Balanced sensitivity/specificity |
| **Viral episomes** | 1.5 | 5 | 100 | Detect small viral circles |
| **Large eccDNA/amplicons** | 3.0 | 15 | 10000 | Focus on large circles |

### Classification Usage

The classification module identifies the genomic composition of detected regions:

```bash
# Basic classification
circontrack-classify circdna.bed reference.fa aligned.bam -o results

# With custom viral contig patterns (for RefSeq viruses)
circontrack-classify circdna.bed reference.fa aligned.bam \
  --viral-patterns NC_ NR_ viral \
  -o results

# Outputs:
# - results_classified.bed (annotated regions)
# - results_summary.txt (statistics)
```

## Output Format

### Detection Output (BED)

```
#chrom  start    end      name        score  strand  method    length  fold_enrichment  confidence
chr19   1000000  1005000  circDNA_1   850    .       coverage  5000    3.5             0.85
NC_029549 0      15234    circDNA_2   950    .       junction  15234   8.2             0.95
```

### Classification Output

```
#chrom  start    end      name        score  strand  type              confidence  reads
chr19   1000000  1005000  circDNA_1   850    .       host_eccDNA       0.85        120
NC_029549 0      15234    circDNA_2   950    .       viral_episome     0.95        89
chr3    5000000  5002000  circDNA_3   900    .       integration_site  0.90        95
```

## Complete Workflow Example

```bash
# 1. Prepare reference (host + viral)
cat mm10.fa viral_refseq.fa > combined_ref.fa
samtools faidx combined_ref.fa

# 2. Align ONT reads
minimap2 -ax map-ont -t 8 combined_ref.fa reads.fastq.gz | \
  samtools sort -@ 8 -o sample.bam
samtools index sample.bam

# 3. Detect circular DNA
circontrack sample.bam combined_ref.fa \
  --min-fold-enrichment 2.0 \
  --min-coverage 10 \
  --min-length 200 \
  -o circdna_detected.bed

# 4. Classify regions (if viral sequences present)
circontrack-classify circdna_detected.bed combined_ref.fa sample.bam \
  --viral-patterns NC_ NR_ \
  -o final_results

# 5. View results
head final_results_classified.bed
cat final_results_summary.txt
```

## Interpreting Results

### Region Types

| Type | Description | Biological Significance |
|------|-------------|------------------------|
| **host_eccDNA** | Extrachromosomal circular DNA from host genome | Normal or cancer-related eccDNA |
| **viral_episome** | Free circular viral DNA | Active viral replication |
| **integration_site** | Viral DNA integrated into host | Potential oncogenic driver |
| **chimeric** | Complex host-viral rearrangement | Genome instability |

### Confidence Scores

- **0.9-1.0**: Very high confidence
- **0.7-0.9**: High confidence  
- **0.5-0.7**: Moderate confidence
- **0.3-0.5**: Low confidence
- **<0.3**: Filtered by default

## Performance Tips

```bash
# For large datasets, process chromosomes in parallel
parallel -j 8 "circontrack sample.bam ref.fa -c chr{} -o chr{}.bed" ::: {1..22} X Y

# Merge results
cat chr*.bed | grep -v "^#" | sort -k1,1 -k2,2n > all_circles.bed

# For high coverage samples, increase thresholds
circontrack sample.bam ref.fa \
  --min-fold-enrichment 5.0 \
  --min-coverage 30 \
  -o high_coverage_circles.bed
```

## Requirements

- Python ≥ 3.7
- pysam ≥ 0.19.0
- numpy ≥ 1.19.0
- scipy ≥ 1.6.0

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No circles detected | Lower `--min-fold-enrichment` to 1.5 |
| Too many false positives | Increase `--min-coverage` to 15-20 |
| Missing small circles | Decrease `--min-length` to 100 |
| Missing viral contigs | Check reference includes viral sequences |
| Classification not working | Ensure viral contigs start with NC_, NR_, or use `--viral-patterns` |

## Citation

If you use CircONTrack in your research, please cite:

```
CircONTrack: ONT-optimized circular DNA detection and classification
[Citation details to be added upon publication]
```

## Support

- GitHub: [https://github.com/samarth8392/CircONTrack](https://github.com/samarth8392/CircONTrack)
- Issues: [https://github.com/samarth8392/CircONTrack/issues](https://github.com/samarth8392/CircONTrack/issues)

## License

MIT License - see [LICENSE](https://github.com/samarth8392/CircONTrack/blob/main/LICENSE) file for details.