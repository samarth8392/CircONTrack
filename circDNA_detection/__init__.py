#!/usr/bin/env python3
"""
CircONTrack: ONT-optimized Circular DNA Detection Pipeline
"""

__version__ = "0.7.2"
__author__ = "Samarth Mathur, PhD"
__email__ = "samarth8392@gmail.com"

# Import main components
from .circular_dna_detector import CircularDNADetector
from .coverage_analyzer import CoverageAnalyzer
from .junction_detector import JunctionDetector
from .split_read_analyzer import SplitReadAnalyzer
from .confidence_scorer import ConfidenceScorer, MultiMethodIntegrator
from .utils import CircularCandidate, filter_candidates_by_confidence, calculate_gc_content

__all__ = [
    'CircularDNADetector',
    'CoverageAnalyzer',
    'JunctionDetector',
    'SplitReadAnalyzer',
    'ConfidenceScorer',
    'MultiMethodIntegrator',
    'CircularCandidate',
    'filter_candidates_by_confidence',
    'calculate_gc_content',
]