## Summary

The CircONTrack classification module (`circontrack-classify`) provides a streamlined, single-BAM approach to characterizing circular DNA regions by:

1. Leveraging long-read advantages to detect chimeric molecules directly
2. Using SA tags and soft clips to identify integration junctions
3. Classifying regions based on read composition patterns
4. Providing confidence scores for each classification

## Output Formats

### Classified BED Format

```
chr19  1000  5000  circDNA_1  850  .  integration_site  0.85  100  45  20  30  5  High Plains wheat mosaic virus|CMV  NC_029549,NC_029550
```

### Summary Report Structure

```
Viral species detected:
  High Plains wheat mosaic virus: 12 regions
  Cytomegalovirus: 8 regions
  Human papillomavirus type 16: 3 regions

Viral contigs (RefSeq accessions):
  NC_029549: 12 regions (High Plains wheat mosaic virus)
  NC_001699: 8 regions (Cytomegalovirus)
  ...
```

## Usage Examples

### Basic Classification
Basic usage with RefSeq viruses:

```bash
bashcircontrack-classify circdna.bed combined_ref.fa combined.bam -o results
```


### Custom viral patterns:
```bash
# If your viruses have different naming patterns
circontrack-classify circdna.bed combined_ref.fa combined.bam \
    --viral-patterns NC_ NR_ phage_ virus_ -o results
```

### Without parsing species names:
```bash
# Just use contig IDs
circontrack-classify circdna.bed combined_ref.fa combined.bam \
    --no-parse-names -o results
```

## Biological Interpretation Guide

| Classification | Biological Meaning | Key Features | Clinical Relevance |
|---------------|-------------------|--------------|-------------------|
| `host_eccDNA` | Normal extrachromosomal circular DNA | >80% host reads, no viral | Background eccDNA |
| `viral_episome` | Free circular viral DNA | >90% viral reads, circular signatures | Active viral replication |
| `viral_dominant` | Mostly viral with some host | 70-90% viral | Possible contamination or early integration |
| `integration_site` | Viral DNA integrated into host genome | Junction reads, chimeric signatures | Potential oncogenic |
| `chimeric` | Complex host-viral rearrangement | Multiple transitions | Genome instability |
| `mixed` | Unclear classification | No dominant type | Requires manual review |

## Validation and Quality Control

### Expected Patterns by Sample Type

| Sample Type | Expected Distribution |
|------------|----------------------|
| Normal tissue | >95% host_eccDNA, <5% other |
| Viral infection | Mix of host_eccDNA and viral_episome |
| Cancer with viral integration | Presence of integration_site regions |
| Cell line with viral vectors | High proportion of viral_episome |

### Quality Metrics

1. **Read depth**: Minimum 10 reads per region for confident classification
2. **Junction support**: ≥2 junction reads for integration site confirmation  
3. **Consistency**: Multiple reads showing same pattern increases confidence

## Troubleshooting

| Issue | Possible Cause | Solution |
|-------|---------------|----------|
| No viral detected | Wrong viral prefix | Check reference contig names |
| All regions "mixed" | Low coverage | Increase sequencing depth |
| False integration sites | Short read contamination | Increase min_mapq threshold |
| Missing episomes | Circular reads not detected | Check aligner settings for supplementary alignments |

