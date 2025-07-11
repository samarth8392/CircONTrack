#!/usr/bin/env python3
"""
Basic import and functionality tests for CircONTrack
"""

import pytest
from circDNA_detection.utils import CircularCandidate, filter_candidates_by_confidence, calculate_gc_content
from circDNA_detection.confidence_scorer import ConfidenceScorer, MultiMethodIntegrator


def test_circular_candidate_creation():
    """Test CircularCandidate dataclass creation"""
    candidate = CircularCandidate(
        chromosome="chr1",
        start=1000,
        end=2000,
        length=0  # Should be calculated automatically
    )
    assert candidate.length == 1000
    assert candidate.confidence_score == 0.0
    assert candidate.detection_method == 'unknown'


def test_gc_content_calculation():
    """Test GC content calculation"""
    assert calculate_gc_content("ATCG") == 0.5
    assert calculate_gc_content("AAAA") == 0.0
    assert calculate_gc_content("GGCC") == 1.0
    assert calculate_gc_content("ATCGN") == 0.5  # N should be excluded
    assert calculate_gc_content("") == 0.0


def test_confidence_filtering():
    """Test filtering candidates by confidence"""
    candidates = [
        CircularCandidate("chr1", 1000, 2000, 1000, confidence_score=0.8),
        CircularCandidate("chr1", 3000, 4000, 1000, confidence_score=0.2),
        CircularCandidate("chr1", 5000, 6000, 1000, confidence_score=0.5),
    ]
    
    filtered = filter_candidates_by_confidence(candidates, min_confidence=0.4)
    assert len(filtered) == 2
    assert all(c.confidence_score >= 0.4 for c in filtered)


def test_confidence_scorer():
    """Test ConfidenceScorer functionality"""
    scorer = ConfidenceScorer()
    
    # Test coverage-based candidate
    coverage_candidate = CircularCandidate(
        chromosome="chr1",
        start=1000,
        end=2000,
        length=1000,
        detection_method="coverage",
        fold_enrichment=3.0,
        coverage_uniformity=0.8,
        mean_coverage=15.0
    )
    
    score = scorer.calculate_confidence(coverage_candidate)
    assert 0.0 <= score <= 1.0
    assert score > 0.5  # Should be reasonably high with these values


def test_multi_method_integrator():
    """Test MultiMethodIntegrator"""
    integrator = MultiMethodIntegrator()
    
    coverage_cands = [
        CircularCandidate("chr1", 1000, 2000, 1000, detection_method="coverage")
    ]
    junction_cands = [
        CircularCandidate("chr1", 1500, 2500, 1000, detection_method="junction")
    ]
    split_cands = []
    
    integrated = integrator.integrate_candidates(coverage_cands, junction_cands, split_cands)
    
    assert len(integrated) == 1  # Should merge overlapping candidates
    assert "coverage" in integrated[0].detection_method
    assert "junction" in integrated[0].detection_method


if __name__ == "__main__":
    pytest.main([__file__, "-v"])