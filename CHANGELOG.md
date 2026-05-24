# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

### Added
- Added optional `circontrack --report` and `--plot-dir` outputs for Markdown candidate reports and candidate summary plots.
- Added `circDNA_detection.reporting` helpers for structured candidate reports, confidence summaries, evidence-support plots, confidence-breakdown plots, and coverage-profile rendering.
- Added `requirements-dev.txt` with pytest and plotting/reporting dependencies for local validation.

### Changed
- Changed candidate plotting/reporting to render from final `CircularCandidate` objects without recomputing detection evidence.
- Changed `circontrack` BAM validation to require an indexed BAM and report missing indexes as invalid input because random-access BAM fetches are used.
- Changed peak-analysis plotting in `coverage_peakout.py` and `coverage_peaks.py` to use deterministic, headless Matplotlib rendering with clearer labels and colorblind-aware colors.

### Fixed
- Fixed coverage flanking-background calculation so the right flank excludes the candidate peak interval; coverage fold-enrichment values may differ for affected coverage candidates.
- Fixed `coverage_peakout.py` and `circontrack_validate.py` peak-file loading so the first non-comment data row is preserved when files use a commented `#chr` header.
- Fixed an invalid escape warning in the `circontrack` CLI banner.

### Documentation
- Added repository architecture documentation covering CLI entry points, package structure, primary data flow, candidate integration, and output generation.
- Added scientific methodology documentation with implemented formulas, LaTeX notation, variable definitions, coordinate conventions, edge cases, limitations, synthetic test-data design, and ASCII pipeline/evidence diagrams.
- Updated README and docs index language to describe implemented candidate-detection behavior and avoid unsupported accuracy claims.
- Added documentation for optional report and plot outputs and clarified the primary BED output schema.

### Tests
- Added deterministic synthetic FASTA/BAM pytest fixtures using pysam, including empty BAMs, linear reads, elevated coverage regions, SA-tag junction reads, split-read-like evidence, coordinate sorting, and BAM indexing.
- Added unit tests for confidence scoring, confidence filtering, candidate merging, CIGAR/SA parsing, coverage flanking background, and BED output schema.
- Added end-to-end and CLI smoke tests covering empty BAMs, chromosome-specific analysis, contig-boundary clamping, missing BAM index handling, and multi-evidence candidate retention.
- Added plot/report tests confirming Markdown sections and PNG output creation.

## [1.0.0] - 2025-07-10

### Added
- Initial release
- Multi-modal circular DNA detection
- ONT-optimized algorithms
- Coverage pattern analysis
- Junction detection
- Split-read analysis
- Command-line interface
- Python API
