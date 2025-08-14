# CircONTrack Assemble Module

## Overview

The `circontrack-assemble` module extracts reads from detected circular DNA regions and prepares them for *de novo* assembly. This is the final step in the CircONTrack pipeline, enabling reconstruction of complete circular DNA sequences.

## Key Features

- **Targeted read extraction** from circular DNA regions
- **Padding support** for better assembly context
- **Junction-spanning read detection** for circular structures
- **Multi-threaded processing** for large datasets
- **Assembly script generation** for Flye and Canu
- **SLURM/Swarm support** for HPC environments

## Complete CircONTrack Pipeline

```mermaid
flowchart LR
    A[ONT Reads] --> B[Align to Reference]
    B --> C[circontrack detect]
    C --> D[Circular DNA Regions]
    D --> E[circontrack-classify]
    E --> F[Classified Regions]
    F --> G[circontrack-assemble]
    G --> H[Assembled Circles]
    
    style C fill:#e3f2fd
    style E fill:#f3e5f5
    style G fill:#e8f5e9
```

## Usage

### Basic Usage

```bash
# Extract reads from detected regions
circontrack-assemble circdna.bed sample.bam -o assembly_prep

# Use classified BED (includes region types)
circontrack-assemble classified.bed sample.bam --classified -o assembly_prep
```

### Advanced Options

```bash
circontrack-assemble circdna.bed sample.bam \
  --min-reads 100 \        # Minimum reads to attempt assembly
  --min-coverage 20 \       # Minimum coverage depth
  --padding 2000 \          # Include 2kb padding around regions
  --assembler flye \        # Use Flye assembler
  --slurm \                 # Create SLURM submission script
  -t 8 \                    # Use 8 threads
  -o assembly_prep
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-reads` | 50 | Minimum reads required for assembly |
| `--min-coverage` | 10 | Minimum coverage depth for assembly |
| `--padding` | 1000 | Extra bases around regions (improves assembly) |
| `--threads` | 4 | Number of parallel threads |
| `--assembler` | flye | Assembler to use (flye or canu) |
| `--slurm` | False | Create SLURM/Swarm submission script |
| `--classified` | False | Input is classified BED with region types |

## Output Structure

```
assembly_prep/
├── reads/                    # Extracted FASTQ files
│   ├── chr1_1000_5000_circDNA_1.fastq
│   ├── chr2_2000_3000_circDNA_2.fastq
│   └── ...
├── info/                     # Detailed region information
│   ├── chr1_1000_5000_circDNA_1.info
│   └── ...
├── extraction_summary.txt    # Overall summary
├── extraction_stats.tsv      # Detailed statistics
├── assembly_jobs.txt         # List of assembly jobs
└── run_flye_assembly.sh     # Assembly script
```

## Region Information Files

Each `.info` file contains:
- Region coordinates and padding
- Read statistics (total, primary, supplementary)
- Junction-spanning read counts
- Coverage depth and N50
- Assembly recommendations

Example:
```
CircONTrack Region Information
==================================================
Region: chr19:1000000-1005000
Region ID: chr19_1000000_1005000_circDNA_1
Region Type: host_eccDNA
Length: 5,000 bp
Padded region: chr19:999000-1006000

Read Statistics:
Reads extracted: 250
Junction-spanning reads: 12
Average read length: 8,500 bp
Estimated coverage: 42.5x

Assembly Recommendation:
✓ Moderate coverage - suitable for assembly
✓ Junction-spanning reads detected - circular structure likely
```

## Assembly Script Generation

### Flye Assembly (Recommended for Circular DNA)

```bash
# Generated script runs:
flye --nano-raw reads/region_1.fastq \
  --out-dir assembly_region_1 \
  --threads 4 \
  --meta \                    # Metagenomic mode for uneven coverage
  --min-overlap 1000 \        # Smaller overlap for circular DNA
  --iterations 2              # Polish twice
```

### For HPC/SLURM Environments

```bash
# Create SLURM submission script
circontrack-assemble circdna.bed sample.bam --slurm -o assembly_prep

# Submit jobs
bash assembly_prep/submit_flye.sh
```

## Complete Workflow Example

```bash
# 1. Detect circular DNA
circontrack sample.bam reference.fa \
  --min-fold-enrichment 2.0 \
  -o circdna.bed

# 2. Classify regions (if viral sequences in reference)
circontrack-classify circdna.bed reference.fa sample.bam \
  -o classified

# 3. Extract reads for assembly
circontrack-assemble classified_classified.bed sample.bam \
  --classified \
  --min-coverage 15 \
  -o assembly_prep

# 4. Run assembly
bash assembly_prep/run_flye_assembly.sh

# 5. Check assembled sequences
for dir in assembly_prep/assembly_*/; do
  echo "Region: $(basename $dir)"
  if [ -f "$dir/assembly.fasta" ]; then
    grep ">" "$dir/assembly.fasta"
    seqkit stats "$dir/assembly.fasta"
  fi
done
```

## Filtering Recommendations

### By Region Type (for classified BED)

| Region Type | Min Reads | Min Coverage | Notes |
|-------------|-----------|--------------|-------|
| host_eccDNA | 50 | 10x | Standard thresholds |
| viral_episome | 30 | 15x | Often high coverage |
| integration_site | 100 | 20x | Complex, needs more reads |
| chimeric | 100 | 15x | Challenging assembly |

### By Size

| Circle Size | Recommended Coverage | Assembly Notes |
|------------|---------------------|----------------|
| <1 kb | 50x+ | May be too small for long reads |
| 1-10 kb | 30x+ | Ideal for ONT assembly |
| 10-100 kb | 20x+ | Good assembly expected |
| >100 kb | 15x+ | May need parameter tuning |

## Tips for Success

1. **Use adequate padding** - Include 1-2kb padding to capture reads extending beyond the circular DNA boundaries

2. **Check junction reads** - Regions with junction-spanning reads are more likely true circles

3. **Prioritize high-coverage regions** - Start with regions having >30x coverage for best results

4. **Viral episomes** - These often have very high coverage and assemble well

5. **Integration sites** - May be challenging; consider manual curation

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No reads extracted | Check BAM file path and index |
| Low coverage for all regions | Verify correct BAM file used for detection |
| Assembly fails | Increase --min-reads threshold |
| Fragmented assemblies | Increase --min-coverage requirement |
| Out of memory | Reduce --threads or use SLURM with more memory |

## Integration with Other Tools

### Verify Circular Structure

After assembly, verify circular structures:

```bash
# Check for circular overlaps
for fasta in assembly_prep/assembly_*/assembly.fasta; do
  nucmer --maxmatch -c 100 -p self $fasta $fasta
  show-coords -r -c -l self.delta | grep -v "^=" | head
done
```

### Polish Assemblies

Polish with original reads:

```bash
# Using Medaka (ONT-specific)
medaka_consensus -i reads.fastq -d assembly.fasta -o polished
```

## Citation

If you use CircONTrack Assemble in your work:

```
CircONTrack: ONT-optimized circular DNA detection and assembly
[Citation to be added]
```