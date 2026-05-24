"""Synthetic pysam BAM tests for CircONTrack detection behavior."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pysam")

from circDNA_detection.circular_dna_detector import CircularDNADetector
from circDNA_detection.coverage_analyzer import CoverageAnalyzer
from circDNA_detection.junction_detector import JunctionDetector
from circDNA_detection.split_read_analyzer import SplitReadAnalyzer
from circDNA_detection.utils import CircularCandidate

from conftest import (
    build_header,
    make_read,
    make_split_junction_read,
    write_bam,
    write_reference,
)


def test_empty_bam_returns_no_candidates(tmp_path, default_contigs):
    reference = write_reference(tmp_path, default_contigs)
    header = build_header(default_contigs)
    bam = write_bam(tmp_path, header, [], name="empty")

    detector = CircularDNADetector(verbose=False, log_level="ERROR")
    output = tmp_path / "empty.bed"
    candidates = detector.detect_circular_dna(str(bam), str(reference), str(output))

    assert candidates == []
    assert output.exists()
    assert output.read_text().startswith("# CircDNA Detection Results")


def test_linear_only_reads_do_not_create_candidates(tmp_path, default_contigs):
    reference = write_reference(tmp_path, default_contigs)
    header = build_header(default_contigs)
    reads = [
        make_read(header, f"linear_{i}", "chr1", 1000 + i * 200, cigar=((0, 120),))
        for i in range(20)
    ]
    bam = write_bam(tmp_path, header, reads, name="linear")

    detector = CircularDNADetector(verbose=False, log_level="ERROR")
    candidates = detector.detect_circular_dna(str(bam), str(reference), chromosome="chr1")

    assert candidates == []


def test_elevated_coverage_region_detected(tmp_path, default_contigs):
    header = build_header(default_contigs)
    reads = []
    window_size = 500

    for window_index, start in enumerate(range(0, default_contigs["chr1"], window_size)):
        read_count = 10 if 5000 <= start < 7000 else 1 + (window_index % 2)
        for j in range(read_count):
            reads.append(
                make_read(
                    header,
                    f"cov_{window_index}_{j}",
                    "chr1",
                    start + 20 + j,
                    cigar=((0, 100),),
                )
            )

    bam = write_bam(tmp_path, header, reads, name="coverage_peak")
    analyzer = CoverageAnalyzer(
        window_sizes=[window_size],
        min_fold_enrichment=1.5,
        min_coverage=5,
        verbose=False,
    )
    candidates = analyzer.detect_coverage_patterns(str(bam), chromosome="chr1")

    assert any(c.start <= 6000 < c.end for c in candidates)
    assert all(c.detection_method == "coverage" for c in candidates)


def test_junction_supporting_reads_detected(tmp_path, default_contigs):
    header = build_header(default_contigs)
    reads = [
        make_split_junction_read(header, f"junction_{i}", "chr1", 100 + i, 250 + i)
        for i in range(3)
    ]
    bam = write_bam(tmp_path, header, reads, name="junction")

    detector = JunctionDetector(min_support=1, min_alignment_length=50, verbose=False)
    candidates = detector.detect_junctions(str(bam), chromosome="chr1")

    assert candidates
    assert candidates[0].junction_support >= 1
    assert candidates[0].start == 0
    assert candidates[0].detection_method == "junction"


def test_split_read_like_evidence_detected(tmp_path, default_contigs):
    header = build_header(default_contigs)
    reads = [
        make_split_junction_read(header, f"split_{i}", "chr1", 100 + i, 250 + i)
        for i in range(3)
    ]
    bam = write_bam(tmp_path, header, reads, name="split")

    analyzer = SplitReadAnalyzer(
        min_support=1,
        min_split_length=50,
        min_split_mapq=20,
        max_distance=1000,
        verbose=False,
    )
    candidates = analyzer.analyze_split_reads(str(bam), chromosome="chr1")

    assert candidates
    assert candidates[0].split_support >= 1
    assert candidates[0].start == 0
    assert candidates[0].detection_method == "split_read"


def test_chromosome_specific_analysis(tmp_path, default_contigs):
    header = build_header(default_contigs)
    reads = [
        make_split_junction_read(header, "chr1_signal", "chr1", 100, 250),
        make_read(header, "chr2_linear", "chr2", 100, cigar=((0, 150),)),
    ]
    bam = write_bam(tmp_path, header, reads, name="chromosome_specific")
    detector = JunctionDetector(min_support=1, min_alignment_length=50, verbose=False)

    assert detector.detect_junctions(str(bam), chromosome="chr1")
    assert detector.detect_junctions(str(bam), chromosome="chr2") == []


def test_contig_boundary_split_candidate_is_clamped(tmp_path, default_contigs):
    header = build_header(default_contigs)
    reads = [
        make_split_junction_read(
            header,
            "boundary_split",
            "chr1",
            primary_start=5,
            supplementary_start=155,
        )
    ]
    bam = write_bam(tmp_path, header, reads, name="boundary")
    analyzer = SplitReadAnalyzer(min_support=1, min_split_length=50, verbose=False)
    candidates = analyzer.analyze_split_reads(str(bam), chromosome="chr1")

    assert candidates
    assert min(c.start for c in candidates) == 0


def test_missing_bam_index_raises_clear_error(tmp_path, default_contigs):
    reference = write_reference(tmp_path, default_contigs)
    header = build_header(default_contigs)
    bam = write_bam(tmp_path, header, [], name="unindexed", index=False)
    detector = CircularDNADetector(verbose=False, log_level="ERROR")

    with pytest.raises(ValueError, match="indexed"):
        detector._validate_inputs(str(bam), str(reference))


def test_invalid_input_handling(tmp_path, default_contigs):
    reference = write_reference(tmp_path, default_contigs)
    detector = CircularDNADetector(verbose=False, log_level="ERROR")

    with pytest.raises(ValueError, match="Invalid BAM file"):
        detector._validate_inputs(str(tmp_path / "missing.bam"), str(reference))


def test_bed_output_schema_correctness(tmp_path):
    output = tmp_path / "results.bed"
    detector = CircularDNADetector(verbose=False, log_level="ERROR")
    candidate = CircularCandidate(
        "chr1",
        10,
        210,
        200,
        confidence_score=0.75,
        detection_method="junction+split_read",
        gc_content=0.5,
    )

    detector._write_output([candidate], str(output))
    lines = output.read_text().splitlines()

    assert lines[1] == "# chr\tstart\tend\tname\tconfidence\tstrand\tmethod\tlength\tgc_content"
    fields = lines[2].split("\t")
    assert fields == [
        "chr1",
        "10",
        "210",
        "circDNA_1",
        "0.750",
        ".",
        "junction+split_read",
        "200",
        "0.500",
    ]


def test_end_to_end_multi_evidence_candidate_retained(tmp_path, default_contigs):
    reference = write_reference(tmp_path, default_contigs)
    header = build_header(default_contigs)
    reads = [
        make_split_junction_read(header, f"multi_{i}", "chr1", 100 + i, 250 + i)
        for i in range(3)
    ]
    bam = write_bam(tmp_path, header, reads, name="multi_evidence")
    output = tmp_path / "multi.bed"

    detector = CircularDNADetector(
        min_confidence=0.3,
        verbose=False,
        log_level="ERROR",
    )
    candidates = detector.detect_circular_dna(
        str(bam),
        str(reference),
        str(output),
        chromosome="chr1",
    )

    assert candidates
    assert "junction" in candidates[0].detection_method
    assert "split_read" in candidates[0].detection_method
    assert candidates[0].confidence_score >= 0.3
    assert len([line for line in output.read_text().splitlines() if not line.startswith("#")]) >= 1


def test_cli_smoke_empty_bam(tmp_path, default_contigs):
    reference = write_reference(tmp_path, default_contigs)
    header = build_header(default_contigs)
    bam = write_bam(tmp_path, header, [], name="cli_empty")
    output = tmp_path / "cli_empty.bed"

    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "circDNA_detection.circular_dna_detector",
            str(bam),
            str(reference),
            "-o",
            str(output),
            "-q",
            "--log-level",
            "ERROR",
            "--min-confidence",
            "0.99",
        ],
        cwd=str(Path(__file__).resolve().parents[1]),
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert output.exists()
    assert "Found 0 high-confidence circular DNA candidates" in result.stdout
