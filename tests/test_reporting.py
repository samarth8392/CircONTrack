"""Tests for report generation, plotting, and peak-output parsing."""

import pytest

pytest.importorskip("pysam")

from circDNA_detection.reporting import (
    plot_coverage_profile,
    write_markdown_report,
    write_plot_bundle,
)
from circDNA_detection.utils import CircularCandidate

from conftest import build_header, make_read, write_bam


def _candidate():
    return CircularCandidate(
        "chr1",
        100,
        500,
        400,
        confidence_score=0.82,
        detection_method="coverage+junction+split_read",
        mean_coverage=18.0,
        fold_enrichment=4.0,
        coverage_uniformity=0.85,
        junction_support=5,
        split_support=4,
        gc_content=0.5,
    )


def test_markdown_report_contains_expected_sections(tmp_path):
    report = tmp_path / "candidate_report.md"
    write_markdown_report([_candidate()], str(report), run_metadata={"chromosome": "chr1"})
    text = report.read_text()

    assert "# CircONTrack Candidate Report" in text
    assert "## Candidate Summary" in text
    assert "## Evidence Table" in text
    assert "## Confidence Score Breakdown" in text
    assert "0-based, half-open" in text


def test_candidate_plot_bundle_creates_files(tmp_path):
    pytest.importorskip("matplotlib")

    paths = write_plot_bundle([_candidate()], str(tmp_path / "plots"))

    assert {path.name for path in paths} == {
        "candidate_confidence.png",
        "candidate_evidence_support.png",
        "confidence_breakdown.png",
    }
    assert all(path.exists() and path.stat().st_size > 0 for path in paths)


def test_coverage_profile_plot_created_from_structured_windows(tmp_path):
    pytest.importorskip("matplotlib")

    windows = [
        {"start": 0, "end": 100, "coverage": 1.0},
        {"start": 100, "end": 200, "coverage": 2.0},
        {"start": 200, "end": 300, "coverage": 8.0},
    ]
    output = tmp_path / "coverage_profile.png"

    plot_coverage_profile(windows, [_candidate()], "chr1", str(output))

    assert output.exists()
    assert output.stat().st_size > 0


def test_peakout_loader_preserves_first_data_row(tmp_path):
    pytest.importorskip("pandas")
    pytest.importorskip("matplotlib")

    from circDNA_detection.coverage_peakout import CircONTrackPeakAnalyzer

    peak_file = tmp_path / "peaks.bed"
    peak_file.write_text(
        "#chr\tstart\tend\tname\tscore\tstrand\tcoverage\tfold_change\tpvalue\tadjusted_pvalue\tread_count\n"
        "chr1\t10\t110\tpeak_1\t100\t.\t10\t2.0\t1e-4\t1e-3\t5\n"
        "chr1\t210\t310\tpeak_2\t80\t.\t8\t1.6\t2e-4\t2e-3\t4\n"
    )

    analyzer = CircONTrackPeakAnalyzer(str(peak_file))

    assert len(analyzer.peaks) == 2
    assert list(analyzer.peaks["name"]) == ["peak_1", "peak_2"]


def test_validator_loader_preserves_first_data_row(tmp_path):
    pytest.importorskip("pandas")

    from circDNA_detection.circontrack_validate import CircONTrackValidator

    peak_file = tmp_path / "peaks.bed"
    peak_file.write_text(
        "#chr\tstart\tend\tname\tscore\tstrand\tcoverage\tfold_change\tpvalue\tadjusted_pvalue\tread_count\n"
        "chr1\t10\t110\tpeak_1\t100\t.\t10\t2.0\t1e-4\t1e-3\t5\n"
        "chr1\t210\t310\tpeak_2\t80\t.\t8\t1.6\t2e-4\t2e-3\t4\n"
    )
    header = build_header({"chr1": 1000})
    bam = write_bam(
        tmp_path,
        header,
        [make_read(header, "read_1", "chr1", 20, cigar=((0, 80),))],
        name="validator",
    )

    validator = CircONTrackValidator(str(bam), peak_file=str(peak_file))

    assert len(validator.peaks) == 2
    assert list(validator.peaks["name"]) == ["peak_1", "peak_2"]
