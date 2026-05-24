# Repository Architecture

This map reflects the current implementation.

## Main Entry Points

Console scripts are declared in `setup.py`:

- `circontrack`: primary circular DNA candidate detector, implemented by `circDNA_detection.circular_dna_detector:main`.
- `circontrack-peaks`: statistical coverage peak caller.
- `circontrack-peakout`: downstream analysis of coverage peak BED output.
- `circontrack-validate`: coverage peak validation with read-level evidence.
- `circontrack-classify`: host/viral classification of candidate intervals.
- `circontrack-assemble`: read extraction and assembly-script preparation.
- `circontrack-viral-episomes`: viral episome workflow.

## Primary `circontrack` Flow

```text
CLI arguments
    |
    v
CircularDNADetector
    |
    +--> validate indexed BAM and FASTA
    |
    +--> CoverageAnalyzer.detect_coverage_patterns
    |
    +--> JunctionDetector.detect_junctions
    |
    +--> SplitReadAnalyzer.analyze_split_reads
    |
    v
MultiMethodIntegrator.integrate_candidates
    |
    v
ConfidenceScorer.calculate_confidence
    |
    v
FASTA GC annotation
    |
    v
filter_candidates_by_confidence
    |
    v
BED writer + optional report/plots
```

## Package Structure

| Path | Role |
|---|---|
| `circDNA_detection/circular_dna_detector.py` | Primary CLI and orchestration. |
| `circDNA_detection/utils.py` | `CircularCandidate`, confidence filter, GC fraction utility. |
| `circDNA_detection/coverage_analyzer.py` | Heuristic coverage-pattern candidate detection. |
| `circDNA_detection/junction_detector.py` | SA-tag, soft-clip, and large-gap junction evidence. |
| `circDNA_detection/split_read_analyzer.py` | SA-tag split-read validation and candidate creation. |
| `circDNA_detection/confidence_scorer.py` | Candidate integration and confidence scoring. |
| `circDNA_detection/reporting.py` | Structured Markdown report and candidate summary plots. |
| `circDNA_detection/coverage_peaks.py` | Separate negative-binomial/Poisson coverage peak caller. |
| `circDNA_detection/coverage_peakout.py` | Peak table summary, filtering, plots, and text report. |
| `circDNA_detection/circontrack_validate.py` | Peak validation with junction and artifact metrics. |
| `circDNA_detection/classify.py` | Host/viral classification using contig-name patterns and SA tags. |
| `circDNA_detection/assemble.py` | Read extraction from BED regions and assembly-script creation. |
| `tests/` | Unit, synthetic BAM, CLI smoke, plotting, and report tests. |
| `docs/` | Usage, architecture, methodology, and module documentation. |

## Intermediate Data Structures

The primary detector passes `CircularCandidate` dataclass instances between modules. The dataclass contains:

- Genomic interval: `chromosome`, `start`, `end`, `length`.
- Detection metadata: `detection_method`, `confidence_score`.
- Coverage evidence: `mean_coverage`, `fold_enrichment`, `coverage_uniformity`.
- Junction evidence: `junction_support`.
- Split-read evidence: `split_support`.
- Sequence annotation: `gc_content`.

Coverage-peak utilities use dictionaries and pandas data frames rather than `CircularCandidate`.

## Evidence Integration

```text
coverage candidate: chr1 [1000, 2000)
junction candidate: chr1 [1500, 2500)
split candidate:    chr1 [1600, 2400)
          |
          v
same chromosome and close/overlapping
          |
          v
merged candidate: chr1 [1000, 2500)
methods: coverage+junction+split_read
evidence: max coverage metrics, summed read support
```

## Output Generation

The primary BED writer emits one row per retained candidate:

```text
# chr start end name confidence strand method length gc_content
```

Optional report and plots are generated from the final candidate objects only. They do not recompute candidate detection or alter confidence scores.

