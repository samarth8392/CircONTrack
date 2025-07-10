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
    
    def __post_init__(self):
        """Calculate derived attributes"""
        if not hasattr(self, 'length') or self.length == 0:
            self.length = self.end - self.start

def merge_overlapping_candidates(candidates: List[CircularCandidate], 
                               max_distance: int = 1000) -> List[CircularCandidate]:
    """Merge overlapping candidates from different detection methods"""
    if not candidates:
        return []
    
    # Sort by chromosome and position
    candidates.sort(key=lambda x: (x.chromosome, x.start))
    
    merged = []
    current = candidates[0]
    
    for candidate in candidates[1:]:
        # Check for overlap or proximity
        if (candidate.chromosome == current.chromosome and
            candidate.start <= current.end + max_distance):
            
            # Merge candidates - combine evidence
            merged_candidate = _merge_two_candidates(current, candidate)
            current = merged_candidate
        else:
            merged.append(current)
            current = candidate
    
    merged.append(current)
    return merged

def _merge_two_candidates(c1: CircularCandidate, c2: CircularCandidate) -> CircularCandidate:
    """Merge two overlapping candidates"""
    # Use broader boundaries
    start = min(c1.start, c2.start)
    end = max(c1.end, c2.end)
    
    # Combine detection methods
    methods = set([c1.detection_method, c2.detection_method])
    method_str = '+'.join(sorted(methods))
    
    # Create merged candidate
    merged = CircularCandidate(
        chromosome=c1.chromosome,
        start=start,
        end=end,
        length=end - start,
        detection_method=method_str,
        confidence_score=max(c1.confidence_score, c2.confidence_score)
    )
    
    # Combine attributes
    if c1.mean_coverage is not None or c2.mean_coverage is not None:
        merged.mean_coverage = max(
            c1.mean_coverage or 0, c2.mean_coverage or 0
        )
    
    if c1.fold_enrichment is not None or c2.fold_enrichment is not None:
        merged.fold_enrichment = max(
            c1.fold_enrichment or 0, c2.fold_enrichment or 0
        )
    
    if c1.coverage_uniformity is not None or c2.coverage_uniformity is not None:
        merged.coverage_uniformity = max(
            c1.coverage_uniformity or 0, c2.coverage_uniformity or 0
        )
    
    if c1.junction_support is not None or c2.junction_support is not None:
        merged.junction_support = (c1.junction_support or 0) + (c2.junction_support or 0)
    
    if c1.split_support is not None or c2.split_support is not None:
        merged.split_support = (c1.split_support or 0) + (c2.split_support or 0)
    
    return merged

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of sequence"""
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return gc_count / len(sequence)

def format_candidate_summary(candidate: CircularCandidate) -> str:
    """Format candidate for summary output"""
    summary = f"{candidate.chromosome}:{candidate.start}-{candidate.end}"
    summary += f" ({candidate.length}bp, {candidate.detection_method})"
    
    details = []
    if candidate.fold_enrichment:
        details.append(f"fold={candidate.fold_enrichment:.1f}")
    if candidate.junction_support:
        details.append(f"junctions={candidate.junction_support}")
    if candidate.split_support:
        details.append(f"splits={candidate.split_support}")
    
    if details:
        summary += f" [{', '.join(details)}]"
    
    return summary

def filter_candidates_by_confidence(candidates: List[CircularCandidate], 
                                   min_confidence: float = 0.3) -> List[CircularCandidate]:
    """Filter candidates by confidence threshold"""
    return [c for c in candidates if c.confidence_score >= min_confidence]

def get_top_candidates(candidates: List[CircularCandidate], 
                      n: int = 10) -> List[CircularCandidate]:
    """Get top N candidates by confidence score"""
    sorted_candidates = sorted(candidates, key=lambda x: x.confidence_score, reverse=True)
    return sorted_candidates[:n]