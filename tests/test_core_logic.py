"""Unit tests for coordinate, scoring, merging, and parsing logic."""

import numpy as np
import pytest

pytest.importorskip("pysam")

from circDNA_detection.confidence_scorer import ConfidenceScorer, MultiMethodIntegrator
from circDNA_detection.coverage_analyzer import CoverageAnalyzer
from circDNA_detection.junction_detector import JunctionDetector
from circDNA_detection.split_read_analyzer import SplitReadAnalyzer
from circDNA_detection.utils import CircularCandidate, filter_candidates_by_confidence


def test_coverage_confidence_formula_matches_implementation():
    candidate = CircularCandidate(
        "chr1",
        100,
        1100,
        1000,
        detection_method="coverage",
        fold_enrichment=3.0,
        coverage_uniformity=0.8,
        mean_coverage=15.0,
    )
    score = ConfidenceScorer().calculate_confidence(candidate)
    expected = (3.0 / 5.0) * 0.4 + 0.8 * 0.3 + (15.0 / 20.0) * 0.3
    assert score == pytest.approx(expected)


def test_low_confidence_candidate_filtered_out():
    candidates = [
        CircularCandidate("chr1", 100, 300, 200, confidence_score=0.29),
        CircularCandidate("chr1", 500, 900, 400, confidence_score=0.30),
    ]
    retained = filter_candidates_by_confidence(candidates, min_confidence=0.30)
    assert retained == [candidates[1]]


def test_high_confidence_multi_evidence_candidate_retained():
    coverage = CircularCandidate(
        "chr1",
        1000,
        2000,
        1000,
        detection_method="coverage",
        fold_enrichment=5.0,
        coverage_uniformity=0.9,
        mean_coverage=25,
    )
    junction = CircularCandidate(
        "chr1",
        1100,
        2100,
        1000,
        detection_method="junction",
        junction_support=8,
    )
    split = CircularCandidate(
        "chr1",
        1200,
        2200,
        1000,
        detection_method="split_read",
        split_support=7,
    )
    integrated = MultiMethodIntegrator(max_overlap_distance=100).integrate_candidates(
        [coverage], [junction], [split]
    )
    assert len(integrated) == 1
    integrated[0].confidence_score = ConfidenceScorer().calculate_confidence(integrated[0])
    retained = filter_candidates_by_confidence(integrated, min_confidence=0.7)
    assert len(retained) == 1
    assert retained[0].detection_method == "coverage+junction+split_read"


def test_candidate_merging_is_chromosome_specific_and_deterministic():
    integrator = MultiMethodIntegrator(max_overlap_distance=50)
    chr1 = CircularCandidate("chr1", 100, 500, 400, detection_method="coverage")
    chr2 = CircularCandidate("chr2", 120, 520, 400, detection_method="junction")
    merged = integrator.integrate_candidates([chr1], [chr2], [])
    assert [(c.chromosome, c.detection_method) for c in merged] == [
        ("chr1", "coverage"),
        ("chr2", "junction"),
    ]


def test_flanking_background_excludes_peak_windows():
    analyzer = CoverageAnalyzer(verbose=False)
    coverage = np.array([1, 1, 10, 10, 1, 1], dtype=float)
    background = analyzer._calculate_flanking_background(
        coverage,
        start_idx=2,
        end_idx=3,
        baseline=99,
    )
    assert background == pytest.approx(1.0)


def test_junction_cigar_reference_length_parsing():
    detector = JunctionDetector(min_alignment_length=1)
    assert detector._get_alignment_length("10S20M5I3D4N2H") == 27


def test_split_cigar_stats_parsing():
    analyzer = SplitReadAnalyzer(min_split_length=1)
    stats = analyzer._parse_cigar_stats("10S20M5I3D4N2H")
    assert stats["ref_length"] == 27
    assert stats["aligned_length"] == 25
    assert stats["clipped"] == 12


def test_sa_tag_coordinates_are_converted_to_zero_based():
    detector = JunctionDetector(min_alignment_length=1)
    parsed = detector._parse_sa_tag("chr1,251,+,150M,60,0;")
    assert parsed[0]["start"] == 250
    assert parsed[0]["end"] == 400
