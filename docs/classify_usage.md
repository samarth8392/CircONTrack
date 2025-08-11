## Summary

The CircONTrack classification module (`circontrack-classify`) provides a streamlined, single-BAM approach to characterizing circular DNA regions by:

1. Leveraging long-read advantages to detect chimeric molecules directly
2. Using SA tags and soft clips to identify integration junctions
3. Classifying regions based on read composition patterns
4. Providing confidence scores for each classification

## Usage Examples

### Basic Classification
```bash
# Minimal command
circontrack-classify circdna.bed combined_ref.fa combined.bam -o results

# Output files:
# - results_classified.bed
# - results_summary.txt
```

### With Custom Parameters
```bash
# Custom viral prefix and quality threshold
circontrack-classify circdna.bed combined_ref.fa combined.bam \
    --viral-prefix "virus_" \
    --min-mapq 30 \
    -o high_quality_results
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

