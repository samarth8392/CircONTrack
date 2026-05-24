# CircONTrack Scientific Methodology and Developer Guide

This page documents the behavior implemented in the current CircONTrack codebase. It separates implemented calculations from intended biological interpretation and marks uncertain behavior where the implementation does not make a claim explicit.

## Overview

CircONTrack is a Python package for detecting candidate circular DNA intervals from Oxford Nanopore aligned reads. The primary `circontrack` command combines three evidence streams: windowed coverage enrichment, junction-like read signatures, and split-read signatures. These streams emit `CircularCandidate` objects, which are merged, scored, optionally annotated with GC content from a FASTA reference, filtered by confidence, and written as BED-like output.

The implementation detects computational signals that can be consistent with circular DNA. It does not include a benchmarked accuracy model in this repository, so documentation should not be read as a claim of high sensitivity or specificity.

## Installation

The package metadata is in `setup.py`. The console scripts are installed by:

```bash
pip install -e .
```

For test development, install the development requirements:

```bash
pip install -r requirements-dev.txt
```

## Quick Start

```bash
circontrack reads.sorted.bam reference.fa -o circular_dna_results.bed
```

Optional report and plot outputs from the final candidate table:

```bash
circontrack reads.sorted.bam reference.fa \
  -o circular_dna_results.bed \
  --report circular_dna_report.md \
  --plot-dir circular_dna_plots
```

Analyze one contig:

```bash
circontrack reads.sorted.bam reference.fa -c chr1 -o chr1_candidates.bed
```

## Inputs and Outputs

Inputs:

- Coordinate-sorted, indexed BAM file. The main detector uses `pysam.AlignmentFile.fetch`, so a `.bai` or `.csi` index is required.
- Indexed or indexable FASTA reference for GC-content annotation through `pysam.FastaFile`.
- Optional chromosome/contig name matching the BAM header.

Primary output:

- BED-like candidate table from `CircularDNADetector._write_output`.

Optional outputs:

- Markdown report from `circDNA_detection.reporting.write_markdown_report`.
- PNG plots from `circDNA_detection.reporting.write_plot_bundle`.
- Coverage-peak BED and plots from the separate `circontrack-peaks` command.
- Peak analysis tables and report from `circontrack-peakout`.
- Peak validation tables from `circontrack-validate`.
- Classification and assembly-preparation outputs from `circontrack-classify` and `circontrack-assemble`.

## CLI Reference

Primary command:

```text
circontrack BAM_FILE REFERENCE_FILE
  -o, --output BED
  --report REPORT.md
  --plot-dir DIR
  -c, --chromosome CONTIG
  -q, --quiet
  --log-level DEBUG|INFO|WARNING|ERROR
  --min-fold-enrichment FLOAT
  --min-coverage INT
  --min-length INT
  --max-length INT
  --min-confidence FLOAT
```

Additional console scripts declared in `setup.py`:

- `circontrack-peaks`: statistical coverage peak caller.
- `circontrack-peakout`: summarize and visualize coverage peak output.
- `circontrack-validate`: validate coverage peaks with junction/artifact metrics.
- `circontrack-classify`: classify candidate intervals using host/viral reference naming.
- `circontrack-assemble`: extract reads from candidate intervals for assembly.
- `circontrack-viral-episomes`: viral episome-oriented workflow.

## Example Commands

```bash
circontrack sample.bam hg38.fa -o sample.circdna.bed --report sample.report.md

circontrack sample.bam hg38.fa -c chr7 -o sample.chr7.bed --min-confidence 0.5

circontrack-peaks sample.bam -r hg38.fa --window-size 500 --plot -o sample.peaks.bed

circontrack-peakout sample.peaks.bed -o peak_analysis

circontrack-validate sample.peaks.bed sample.bam -o sample.validated

circontrack-assemble sample.circdna.bed sample.bam -o assembly_prep
```

## Algorithm Overview

The primary detector orchestrates independent evidence modules and integrates their candidates:

```text
BAM + FASTA
    |
    v
[input validation]
    |
    +--> [coverage_analyzer.CoverageAnalyzer]
    |
    +--> [junction_detector.JunctionDetector]
    |
    +--> [split_read_analyzer.SplitReadAnalyzer]
             |
             v
      [confidence_scorer.MultiMethodIntegrator]
             |
             v
      [confidence_scorer.ConfidenceScorer]
             |
             v
      [FASTA GC annotation]
             |
             v
      [confidence filtering]
             |
             v
      BED + optional Markdown report + plots
```

## Pipeline Orchestration

Implemented in `circDNA_detection.circular_dna_detector.CircularDNADetector.detect_circular_dna`:

1. Validate BAM and FASTA inputs.
2. Run coverage pattern analysis.
3. Run junction detection.
4. Run split-read analysis.
5. Merge overlapping/proximal candidates across methods.
6. Recompute confidence scores on integrated candidates.
7. Fetch FASTA sequence for each interval and calculate GC fraction.
8. Retain candidates with `confidence_score >= min_confidence`.
9. Sort retained candidates by descending confidence.
10. Write BED, report, and plots if requested.

## Data-Flow Diagram

```text
                 +------------------+
                 | indexed BAM file |
                 +------------------+
                    |       |      |
                    |       |      |
                    v       v      v
              coverage  junction  split reads
                 |        |          |
                 +--------+----------+
                          |
                          v
               CircularCandidate objects
                          |
                          v
                  overlap/proximity merge
                          |
                          v
               confidence score assignment
                          |
                          v
              FASTA sequence / GC fraction
                          |
                          v
                threshold filtering
                          |
                          v
          BED table + Markdown report + PNG plots
```

## ASCII-Art Evidence Diagrams

Read evidence interpretation:

```text
primary alignment:        [==========]---------- clipped query
SA tag alignment:                 clipped query ----------[==========]
                                      |
                                      v
                         endpoint proximity on same contig
                                      |
                                      v
                              junction/split evidence
```

Coverage candidate detection:

```text
window coverage:  1 2 1 2 9 10 9 2 1
background:       1 1 1 2 2  2 2 2 1
adjusted:         0 1 0 0 7  8 7 0 0
                          ^^^^^
                          peak windows -> candidate interval [s, e)
```

Confidence scoring and filtering:

```text
coverage metrics ----+
junction support ----+--> normalized components --> final confidence S
split-read support --+                                |
                                                      v
                                            keep if S >= threshold
                                                      |
                                                      v
                                                BED/report/plots
```

Output generation:

```text
CircularCandidate list
        |
        +--> BED rows
        |
        +--> Markdown report
        |
        +--> confidence/evidence PNG plots
```

## Core Logic Map

| Area | Module | Main classes/functions | Notes |
|---|---|---|---|
| CLI and orchestration | `circular_dna_detector.py` | `main`, `CircularDNADetector.detect_circular_dna` | Primary `circontrack` entry point. |
| BAM parsing | `coverage_analyzer.py`, `junction_detector.py`, `split_read_analyzer.py` | `pysam.AlignmentFile.fetch` loops | Uses random-access fetches; index required. |
| FASTA handling | `circular_dna_detector.py`, `utils.py` | `_calculate_gc_content`, `calculate_gc_content` | Fetches `[start, end)` sequence from FASTA. |
| Coverage calculation | `coverage_analyzer.py` | `_calculate_windowed_coverage` | Counts reads per fixed window after MAPQ and secondary filters. |
| Candidate region detection | `coverage_analyzer.py` | `_find_coverage_peaks`, `_find_peak_boundaries` | Uses local background subtraction and SciPy `find_peaks`. |
| Junction detection | `junction_detector.py` | `_find_sa_junctions`, `_find_softclip_junctions`, `_find_chimeric_junctions` | Combines SA tags, soft clips, and large CIGAR gaps/skips. |
| Split-read evidence | `split_read_analyzer.py` | `_collect_split_alignments`, `_validate_split_patterns` | SA-tag driven split-read pattern detection. |
| Candidate merging | `confidence_scorer.py` | `MultiMethodIntegrator` | Merges overlapping or nearby intervals across methods. |
| Confidence scoring | `confidence_scorer.py` | `ConfidenceScorer` | Weighted normalized evidence score, bounded at 1.0. |
| Filtering | `utils.py` | `filter_candidates_by_confidence` | Inclusive threshold: `>= min_confidence`. |
| BED writing | `circular_dna_detector.py`, `coverage_peaks.py` | `_write_output`, `CoveragePeakCaller.write_bed` | Writes comment headers and tab-separated rows. |
| Plotting | `reporting.py`, `coverage_peaks.py`, `coverage_peakout.py` | `write_plot_bundle`, `plot_coverage_peaks`, `plot_comprehensive_analysis` | Headless Matplotlib rendering. |
| Report generation | `reporting.py`, `coverage_peakout.py` | `write_markdown_report`, `generate_report` | Candidate report and peak-analysis report. |

## Coverage-Pattern Detection

Implemented in `CoverageAnalyzer`.

Read filters:

- Skip secondary alignments.
- Require `mapping_quality >= 20`.

Window coverage:

\[
C_i = \frac{N_i}{e_i - s_i} W
\]

Variables:

- \(i\): window index.
- \(s_i, e_i\): 0-based half-open window coordinates.
- \(W\): configured `window_size`.
- \(N_i\): number of reads fetched in the window after filters.
- \(C_i\): window-normalized read count. The unit is reads per configured window size, not per-base depth.

Threshold:

\[
T = \max(\tilde{C} + 3 \operatorname{MAD}(C),\ \tilde{C} f_{\min},\ C_{\min})
\]

Variables:

- \(\tilde{C}\): median of all window coverages on the contig.
- \(\operatorname{MAD}(C)\): median absolute deviation.
- \(f_{\min}\): `min_fold_enrichment`, default 1.5 in `circontrack`.
- \(C_{\min}\): `min_coverage`, default 5.

Local background:

\[
B_i = \operatorname{median}\{C_j: i - 10 \le j \le i + 10\}
\]

Adjusted coverage:

\[
A_i = C_i - B_i
\]

SciPy `find_peaks` is called on \(A_i\) with height \(T - \tilde{C}\), distance `max(5, 500 // window_size)`, and width `max(2, 200 // window_size)`.

Candidate interval:

\[
r = [s, e)
\]

where \(s\) is the start coordinate of the left boundary window and \(e\) is the right boundary window start plus `window_size`.

Region metrics:

\[
\bar{C}_r = \operatorname{median}\{C_i: i \in r\}
\]

\[
F_r = \frac{\bar{C}_r}{\max(B_{\mathrm{flank}}, 1)}
\]

\[
U_r = 1 - \frac{\operatorname{MAD}(C_r)}{\max(\bar{C}_r, 1)}
\]

Variables:

- \(F_r\): fold enrichment.
- \(U_r\): coverage uniformity.
- \(B_{\mathrm{flank}}\): median of non-overlapping left and right flanking windows. A previous implementation included the right boundary peak window in the right flank; this is corrected so flanks exclude the candidate interval.

Edge cases:

- Contigs shorter than `window_size * 10` return no coverage candidates.
- If median coverage is below 1 or MAD is zero, coverage detection returns no candidates.
- Candidate length is hard-coded to 200-100000 bp inside coverage detection.

## Junction Detection

Implemented in `JunctionDetector`.

Evidence sources:

```text
read with SA tag
    primary alignment end/start
    supplementary alignment start/end
        |
        v
 endpoint proximity test
        |
        v
 junction position cluster
```

SA tag coordinates are converted from SAM 1-based positions to 0-based positions:

\[
s_{\mathrm{SA}} = p_{\mathrm{SA}} - 1
\]

Reference alignment length from CIGAR:

\[
L_{\mathrm{ref}} = \sum_{\mathrm{op} \in \{M,D,N,=,X\}} l_{\mathrm{op}}
\]

Circular-pattern proximity tests use absolute endpoint distances less than `max_junction_distance`, default 1000 bp:

\[
|e_{\mathrm{primary}} - s_{\mathrm{supp}}| < d_{\max}
\]

or analogous start-to-end and inverted-orientation endpoint comparisons.

Junction clusters are formed by sorted position, adding a new position to the current cluster when:

\[
p_i - \bar{p}_{\mathrm{cluster}} \le d_{\max}
\]

Boundary estimate:

\[
r = [\max(0, \min(P) - 500),\ \max(P) + 500)
\]

where \(P\) contains junction positions and any available primary, supplementary, or span positions.

Edge cases:

- Unique read names determine junction support.
- Candidates require `junction_support >= min_support`, default 3.
- Candidate length is hard-coded to 200-100000 bp.
- Soft-clip-only junction evidence does not realign clipped sequence; it records large clips as positional support.

## Split-Read Detection

Implemented in `SplitReadAnalyzer`.

Supplementary alignments are parsed from SA tags and filtered:

- Same chromosome as primary.
- MAPQ at least `min_split_mapq`, default 20.
- Reference length at least `min_split_length`, default 100 bp.
- If edit distance is present, \(NM / L_{\mathrm{aligned}} \le 0.15\).
- Primary plus supplementary reference lengths cover at least `min_overlap_ratio` of inferred read length, default 0.8.
- Primary/supplementary genomic overlap greater than 50 bp is rejected.

CIGAR statistics:

\[
L_{\mathrm{ref}} = \sum_{\mathrm{op} \in \{M,D,N,=,X\}} l_{\mathrm{op}}
\]

\[
L_{\mathrm{aligned}} = \sum_{\mathrm{op} \in \{M,I,=,X\}} l_{\mathrm{op}}
\]

\[
L_{\mathrm{clip}} = \sum_{\mathrm{op} \in \{S,H\}} l_{\mathrm{op}}
\]

Split quality:

\[
Q_{\mathrm{mapq}} = \frac{(\mathrm{MAPQ}_{p} + \mathrm{MAPQ}_{s}) / 2}{60}
\]

\[
Q_{\mathrm{len}} = \min\left(1,\frac{L_p + L_s}{1000}\right)
\]

\[
Q_{\mathrm{NM}} = 1 - \frac{NM_p / L_p + NM_s / L_s}{2}
\]

\[
Q_{\mathrm{split}} = \frac{Q_{\mathrm{mapq}} + Q_{\mathrm{len}} + Q_{\mathrm{NM}}}{3}
\]

Split candidates are grouped by 500 bp bins:

\[
b = 500 \left\lfloor \frac{p_{\mathrm{junction}}}{500} \right\rfloor
\]

Candidate bounds:

\[
r = [\max(0, \min(P)-200),\ \max(P)+200)
\]

where \(P\) contains primary and supplementary start/end positions.

The split-read confidence score is first calculated by `ConfidenceScorer`, then multiplied by mean split quality:

\[
S_{\mathrm{split, adjusted}} = S_{\mathrm{split}} Q_{\mathrm{split}}
\]

## Candidate Integration

Implemented in `MultiMethodIntegrator`.

Candidates are sorted by chromosome and start coordinate. Two candidates merge when they are on the same chromosome and overlap or lie within `max_overlap_distance`, default 1000 bp:

\[
s_2 \le e_1 + d_{\max} \land s_1 \le e_2 + d_{\max}
\]

Merged bounds:

\[
r_{\mathrm{merged}} = [\min(s_1,s_2),\ \max(e_1,e_2))
\]

Merged evidence:

- Detection methods are de-duplicated and alphabetically joined with `+`.
- Coverage metrics use the maximum available value.
- Junction and split support counts are summed.

## Confidence Scoring

Implemented in `ConfidenceScorer`.

Coverage score:

\[
S_c =
0.4 \min\left(\frac{F_r}{5},1\right)
+ 0.3 \max(U_r,0)
+ 0.3 \min\left(\frac{\bar{C}_r}{20},1\right)
\]

Junction score:

\[
S_j = 0.8 \min\left(\frac{J}{10},1\right) + 0.2
\]

Split-read score:

\[
S_s = 0.8 \min\left(\frac{R}{10},1\right) + 0.2
\]

Variables:

- \(F_r\): fold enrichment.
- \(U_r\): coverage uniformity.
- \(\bar{C}_r\): median coverage in candidate region.
- \(J\): unique junction-supporting reads.
- \(R\): unique split-supporting reads.

Two-method integration:

\[
S = \min\left(\operatorname{mean}(S_m) + 0.2 \min(S_m), 1\right)
\]

Three-method integration:

\[
S = \min\left(\operatorname{mean}(S_m) + 0.3 \min(S_m), 1\right)
\]

where \(S_m\) are the per-method scores for methods in the merged detection-method string.

## Calculation Reference

### Window coverage

Plain language: counts reads overlapping a fixed window and normalizes to the configured window size.

Formula:

\[
C_i = \frac{N_i}{e_i - s_i} W
\]

Implementation: `coverage_analyzer.CoverageAnalyzer._calculate_windowed_coverage`.

Expected input: indexed BAM, chromosome, window size, contig length.

Expected output: list of dictionaries with `start`, `end`, `coverage`, and `raw_count`.

Edge cases: zero-length windows are avoided by window construction; secondary reads are skipped; MAPQ below 20 is skipped.

### Coverage fold enrichment

Plain language: compares median coverage inside a candidate peak with non-overlapping flanking background windows.

Formula:

\[
F_r = \frac{\operatorname{median}(C_r)}{\max(B_{\mathrm{flank}}, 1)}
\]

Implementation: `coverage_analyzer.CoverageAnalyzer._find_coverage_peaks`.

Expected input: per-window coverage array and peak boundaries.

Expected output: candidate `fold_enrichment`.

Edge cases: denominator is floored at 1; missing flanks fall back to baseline; contig-start candidates use only available flanks.

### Coverage uniformity

Plain language: scores whether coverage inside a candidate is relatively even.

Formula:

\[
U_r = 1 - \frac{\operatorname{MAD}(C_r)}{\max(\operatorname{median}(C_r), 1)}
\]

Implementation: `coverage_analyzer.CoverageAnalyzer._find_coverage_peaks`.

Expected input: coverage values inside candidate windows.

Expected output: candidate `coverage_uniformity`.

Edge cases: denominator is floored at 1; values can be negative before confidence scoring clamps the contribution at zero.

### GC content

Plain language: fraction of non-N FASTA bases in the interval that are G or C.

Formula:

\[
G_r = \frac{n_G + n_C}{L - n_N}
\]

Implementation: `utils.calculate_gc_content`, called by `CircularDNADetector._calculate_gc_content`.

Expected input: FASTA sequence string.

Expected output: float between 0 and 1.

Edge cases: empty sequence or all-N sequence returns 0.0.

### Junction support

Plain language: number of unique read names supporting a junction cluster.

Formula:

\[
J = |\{q: q \in \mathrm{reads}(\mathrm{cluster})\}|
\]

Implementation: `junction_detector.JunctionDetector._create_junction_candidates`.

Expected input: clustered junction dictionaries.

Expected output: candidate `junction_support`.

Edge cases: repeated alignments with the same query name count once; support below `min_support` is discarded.

### Split-read support

Plain language: number of unique read names supporting a split-read junction bin.

Formula:

\[
R = |\{q: q \in \mathrm{reads}(\mathrm{bin})\}|
\]

Implementation: `split_read_analyzer.SplitReadAnalyzer._create_split_candidates`.

Expected input: validated split-read dictionaries grouped by junction bin.

Expected output: candidate `split_support`.

Edge cases: bins with fewer than `min_support` split records are discarded; candidate start is clamped to zero.

### Confidence score

Plain language: normalized evidence score used for ranking and threshold filtering.

Formula:

\[
S_c =
0.4 \min\left(\frac{F_r}{5},1\right)
+ 0.3 \max(U_r,0)
+ 0.3 \min\left(\frac{\bar{C}_r}{20},1\right)
\]

\[
S_j = 0.8 \min\left(\frac{J}{10},1\right) + 0.2
\]

\[
S_s = 0.8 \min\left(\frac{R}{10},1\right) + 0.2
\]

Implementation: `confidence_scorer.ConfidenceScorer`.

Expected input: `CircularCandidate` with one or more evidence fields.

Expected output: `confidence_score` in `[0, 1]`.

Edge cases: unknown methods return 0.1; missing evidence components contribute zero except junction/split base confidence when their method is used.

### Coverage peak p-value

Plain language: tests whether a coverage window is enriched relative to a fitted background distribution in the separate `circontrack-peaks` workflow.

Formula for overdispersed windows:

\[
p = \frac{\mu}{\sigma^2}
\]

\[
n = \frac{\mu p}{1-p}
\]

\[
P_i = \Pr(X \ge C_i \mid X \sim \operatorname{NegBin}(n,p))
\]

If variance is not greater than the mean, the implementation uses a Poisson survival function:

\[
P_i = \Pr(X \ge C_i \mid X \sim \operatorname{Poisson}(\mu))
\]

Implementation: `coverage_peaks.CoveragePeakCaller.fit_negative_binomial` and `calculate_pvalues`.

Expected input: window coverage values from `calculate_coverage_windows`.

Expected output: unadjusted p-values and Benjamini-Hochberg adjusted p-values.

Edge cases: fewer than 10 nonzero coverage windows raises an insufficient-data error for distribution fitting; windows with zero background can produce fold change 0 in peak filtering.

## Filtering Thresholds

Final candidate filtering is implemented by `filter_candidates_by_confidence`:

\[
\mathrm{retain}(r) = S_r \ge S_{\min}
\]

Default `S_min` is 0.3 for `circontrack`.

Important implementation note: `CircularDNADetector` accepts `--min-length` and `--max-length`, but the final orchestration layer does not apply an additional length filter. Individual evidence modules use their own hard-coded 200-100000 bp candidate length checks. This is a fragile behavior rather than a documented biological threshold.

## Output BED Schema

Primary `circontrack` output:

| Column | Name | Meaning |
|---:|---|---|
| 1 | `chr` | Candidate chromosome/contig. |
| 2 | `start` | 0-based inclusive start. |
| 3 | `end` | 0-based exclusive end. |
| 4 | `name` | `circDNA_N` rank after confidence sorting. |
| 5 | `confidence` | Final confidence score rounded to three decimals. |
| 6 | `strand` | Always `.` in current implementation. |
| 7 | `method` | Detection method or merged method string. |
| 8 | `length` | `end - start`. |
| 9 | `gc_content` | GC fraction rounded to three decimals, or `NA`. |

Coverage peak output from `circontrack-peaks`:

```text
#chr start end name score strand coverage fold_change pvalue adjusted_pvalue read_count
```

Coordinates are 0-based half-open intervals.

## Plot and Report Interpretation

`circDNA_detection.reporting` renders final candidate summaries from already-computed `CircularCandidate` objects:

- Candidate confidence plot: final confidence score per retained interval.
- Evidence-support plot: junction and split-read support counts per retained interval.
- Confidence-breakdown plot: implemented score components contributing to each candidate.
- Markdown report: run configuration, candidate table, confidence components, and coordinate convention.

The plots visualize the output table. They are not independent validation and should not be interpreted as statistical evidence beyond the encoded metrics.

`coverage_peaks.py` plots coverage windows, highlighted peaks, p-value track, and fold-change track from precomputed windows. Rendering is headless through Matplotlib `Agg`.

## Synthetic Test-Data Design

The pytest suite generates BAM and FASTA inputs at runtime using deterministic helper functions in `tests/conftest.py`.

Synthetic signals:

- Empty BAM: indexed BAM with no reads.
- Linear-only BAM: reads with simple `M` CIGAR operations and no SA tags.
- Coverage enrichment: alternating low background read counts plus a fixed high-count window block.
- Junction-like reads: primary alignment with soft clip and an SA tag whose supplementary segment maps near the primary endpoint.
- Split-read-like reads: the same SA-tag structure, validated through split-read filters.
- Boundary cases: SA-tag signal near contig start, requiring candidate start clamp to 0.
- Missing-index case: coordinate-sorted BAM intentionally left without `.bai`.

No binary BAM files are committed. Test BAMs are coordinate-sorted and indexed programmatically where required.

## Developer Guide

Development priorities:

- Keep scientific behavior explicit and covered by tests before refactoring.
- Prefer small functions that accept structured objects or dictionaries over functions that mix computation, file I/O, and plotting.
- Use deterministic sorting by chromosome/start and candidate confidence.
- Treat BAM coordinates, BED coordinates, and SA-tag coordinates explicitly.
- Avoid adding new runtime dependencies unless the repository already uses them or the behavior cannot be tested otherwise.

Recommended local checks:

```bash
python -m compileall circDNA_detection tests
python -m pytest -q
```

## Testing Guide

The tests cover:

- Empty BAM behavior.
- Linear-only reads.
- Elevated coverage region detection.
- Junction-supporting reads.
- Split-read-like evidence.
- Low-confidence filtering.
- High-confidence multi-evidence retention.
- Chromosome-specific analysis.
- Contig-boundary clamping.
- Missing BAM index error reporting.
- Invalid input handling.
- BED schema correctness.
- CLI smoke behavior.
- Report and plot file creation.
- Peak parser preservation of first data row.

If `pytest` or runtime dependencies are missing, install:

```bash
pip install -r requirements-dev.txt
```

## Limitations and Assumptions

- The repository does not include benchmark data or validated accuracy estimates.
- Coverage detection is sensitive to window size, median/MAD behavior, and local background choice.
- Junction and split-read detection rely on CIGAR and SA-tag heuristics; clipped bases are not re-aligned by this implementation.
- The final detector does not currently apply the CLI `--min-length` and `--max-length` thresholds beyond module-specific hard-coded checks.
- The split-read quality multiplier is applied inside split-read candidate creation and can be recomputed later if candidates are merged, so merged multi-method score behavior should be interpreted as implementation-specific.
- Repetitive regions, low-complexity sequence, mapping artifacts, and reference errors can create signals that resemble circular DNA evidence.

## Troubleshooting

`Invalid BAM file: BAM file must be coordinate-sorted and indexed`:

```bash
samtools sort -o reads.sorted.bam reads.bam
samtools index reads.sorted.bam
```

`Chromosome 'X' not found in BAM header`:

Check contig names:

```bash
samtools idxstats reads.sorted.bam
```

No candidates:

- Lower `--min-confidence` only for exploratory analysis.
- Check whether the BAM has adequate mapped reads on the requested contig.
- Confirm that expected circular junction-supporting reads contain SA tags or informative clipping.

Plotting errors:

- Install Matplotlib in the environment.
- Use a non-interactive/headless environment; CircONTrack plotting uses `Agg`.

Unexpected peak table row counts:

- Coverage peak files should retain all non-comment data rows. The loader explicitly assigns column names because the header line begins with `#`.
