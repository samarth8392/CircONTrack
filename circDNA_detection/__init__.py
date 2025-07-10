#!/usr/bin/env python3
"""
ONT-Optimized Circular DNA Detection Package
Combines coverage patterns, junction detection, and split-read analysis
"""

__version__ = "1.0.0"
__author__ = "Samarth Mathur, PhD"
__email__ = "samarth8392@gmail.com"

from .circular_dna_detector import CircularDNADetector
from .coverage_analyzer import CoverageAnalyzer
from .junction_detector import JunctionDetector
from .split_read_analyzer import SplitReadAnalyzer
from .utils import CircularCandidate, merge_overlapping_candidates, calculate_gc_content

__all__ = [
    'CircularDNADetector',
    'CoverageAnalyzer', 
    'JunctionDetector',
    'SplitReadAnalyzer',
    'CircularCandidate',
    'merge_overlapping_candidates',
    'calculate_gc_content'
]