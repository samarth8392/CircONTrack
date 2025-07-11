#!/usr/bin/env python3
"""
Shared utilities for circular DNA detection
"""

from dataclasses import dataclass
from typing import List, Optional


@dataclass
class CircularCandidate:
    """Unified circular DNA candidate representation"""
    chromosome: str
    start: int
    end: int
    length: int
    confidence_score: float = 0.0
    detection_method: str = 'unknown'
    
    # Coverage-based attributes
    mean_coverage: Optional[float] = None
    fold_enrichment: Optional[float] = None
    coverage_uniformity: Optional[float] = None
    
    # Junction-based attributes
    junction_support: Optional[int] = None
    
    # Split-read attributes
    split_support: Optional[int] = None
    
    # Sequence attributes
    gc_content: Optional[float] = None
    
    def __post_init__(self):
        """Calculate derived attributes"""
        if self.length == 0:
            self.length = self.end - self.start


def filter_candidates_by_confidence(candidates: List[CircularCandidate], 
                                   min_confidence: float = 0.3) -> List[CircularCandidate]:
    """Filter candidates by confidence threshold"""
    return [c for c in candidates if c.confidence_score >= min_confidence]


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of sequence"""
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence) - sequence.count('N')  # Exclude N bases
    
    if total_count == 0:
        return 0.0
    
    return gc_count / total_count