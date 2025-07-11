#!/usr/bin/env python3
"""
Confidence Scoring and Multi-Method Integration for Circular DNA Detection
"""

import numpy as np
from typing import List
from .utils import CircularCandidate


class ConfidenceScorer:
    """
    Centralized confidence scoring system for CircDNA candidates
    """
    
    def __init__(self):
        # Scoring weights for different evidence types
        self.weights = {
            'coverage': {
                'fold_enrichment': 0.4,
                'coverage_uniformity': 0.3,
                'mean_coverage': 0.3
            },
            'junction': {
                'junction_support': 0.8,
                'base_confidence': 0.2
            },
            'split_read': {
                'split_support': 0.8,
                'base_confidence': 0.2
            }
        }
        
        # Reference thresholds for normalization
        self.thresholds = {
            'max_fold_enrichment': 5.0,
            'max_coverage': 20.0,
            'max_junction_support': 10.0,
            'max_split_support': 10.0
        }
    
    def calculate_confidence(self, candidate: CircularCandidate) -> float:
        """
        Calculate confidence score based on available evidence
        Returns normalized score between 0.0 and 1.0
        """
        if candidate.detection_method == 'coverage':
            return self._calculate_coverage_confidence(candidate)
        elif candidate.detection_method == 'junction':
            return self._calculate_junction_confidence(candidate)
        elif candidate.detection_method == 'split_read':
            return self._calculate_split_read_confidence(candidate)
        elif '+' in candidate.detection_method:
            return self._calculate_multi_method_confidence(candidate)
        else:
            return 0.1  # Minimal confidence for unknown methods
    
    def _calculate_coverage_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for coverage-based detection"""
        score = 0.0
        weights = self.weights['coverage']
        
        # Fold enrichment component
        if candidate.fold_enrichment is not None:
            fold_score = min(candidate.fold_enrichment / self.thresholds['max_fold_enrichment'], 1.0)
            score += fold_score * weights['fold_enrichment']
        
        # Coverage uniformity component
        if candidate.coverage_uniformity is not None:
            uniformity_score = max(0, candidate.coverage_uniformity)  # Can be negative
            score += uniformity_score * weights['coverage_uniformity']
        
        # Mean coverage component
        if candidate.mean_coverage is not None:
            coverage_score = min(candidate.mean_coverage / self.thresholds['max_coverage'], 1.0)
            score += coverage_score * weights['mean_coverage']
        
        return min(score, 1.0)
    
    def _calculate_junction_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for junction-based detection"""
        score = 0.0
        weights = self.weights['junction']
        
        # Junction support component
        if candidate.junction_support is not None:
            support_score = min(candidate.junction_support / self.thresholds['max_junction_support'], 1.0)
            score += support_score * weights['junction_support']
        
        # Base confidence for junction detection
        score += weights['base_confidence']
        
        return min(score, 1.0)
    
    def _calculate_split_read_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for split-read detection"""
        score = 0.0
        weights = self.weights['split_read']
        
        # Split support component
        if candidate.split_support is not None:
            support_score = min(candidate.split_support / self.thresholds['max_split_support'], 1.0)
            score += support_score * weights['split_support']
        
        # Base confidence for split-read detection
        score += weights['base_confidence']
        
        return min(score, 1.0)
    
    def _calculate_multi_method_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for multi-method detection"""
        methods = candidate.detection_method.split('+')
        
        # Calculate individual method scores
        method_scores = []
        
        for method in methods:
            # Create temporary candidate with single method for scoring
            temp_candidate = CircularCandidate(
                chromosome=candidate.chromosome,
                start=candidate.start,
                end=candidate.end,
                length=candidate.length,
                detection_method=method,
                mean_coverage=candidate.mean_coverage,
                fold_enrichment=candidate.fold_enrichment,
                coverage_uniformity=candidate.coverage_uniformity,
                junction_support=candidate.junction_support,
                split_support=candidate.split_support
            )
            method_scores.append(self.calculate_confidence(temp_candidate))
        
        # Combine scores with bonus for multi-method support
        if len(method_scores) == 1:
            return method_scores[0]
        elif len(method_scores) == 2:
            # Two methods: average + bonus
            base_score = np.mean(method_scores)
            bonus = 0.2 * min(method_scores)  # Bonus for multi-method support
            return min(base_score + bonus, 1.0)
        else:
            # Three methods: weighted average + larger bonus
            base_score = np.mean(method_scores)
            bonus = 0.3 * min(method_scores)
            return min(base_score + bonus, 1.0)


class MultiMethodIntegrator:
    """
    Integrates candidates from multiple detection methods
    """
    
    def __init__(self, max_overlap_distance=1000):
        self.max_overlap_distance = max_overlap_distance
    
    def integrate_candidates(self, coverage_candidates: List[CircularCandidate],
                           junction_candidates: List[CircularCandidate],
                           split_candidates: List[CircularCandidate]) -> List[CircularCandidate]:
        """
        Integrate candidates from multiple detection methods
        """
        all_candidates = coverage_candidates + junction_candidates + split_candidates
        
        if not all_candidates:
            return []
        
        # Sort by chromosome and position
        all_candidates.sort(key=lambda x: (x.chromosome, x.start))
        
        # Merge overlapping candidates
        merged_candidates = self._merge_overlapping_candidates(all_candidates)
        
        return merged_candidates
    
    def _merge_overlapping_candidates(self, candidates: List[CircularCandidate]) -> List[CircularCandidate]:
        """Merge overlapping candidates from different methods"""
        if not candidates:
            return []
        
        merged = []
        current = candidates[0]
        
        for candidate in candidates[1:]:
            if self._should_merge(current, candidate):
                current = self._merge_two_candidates(current, candidate)
            else:
                merged.append(current)
                current = candidate
        
        merged.append(current)
        return merged
    
    def _should_merge(self, c1: CircularCandidate, c2: CircularCandidate) -> bool:
        """Check if two candidates should be merged"""
        if c1.chromosome != c2.chromosome:
            return False
        
        # Check for overlap or proximity
        overlap = (c2.start <= c1.end + self.max_overlap_distance and
                  c1.start <= c2.end + self.max_overlap_distance)
        
        return overlap
    
    def _merge_two_candidates(self, c1: CircularCandidate, c2: CircularCandidate) -> CircularCandidate:
        """Merge two overlapping candidates"""
        # Use broader boundaries
        start = min(c1.start, c2.start)
        end = max(c1.end, c2.end)
        
        # Combine detection methods
        methods = set()
        if '+' in c1.detection_method:
            methods.update(c1.detection_method.split('+'))
        else:
            methods.add(c1.detection_method)
        
        if '+' in c2.detection_method:
            methods.update(c2.detection_method.split('+'))
        else:
            methods.add(c2.detection_method)
        
        method_str = '+'.join(sorted(methods))
        
        # Create merged candidate
        merged = CircularCandidate(
            chromosome=c1.chromosome,
            start=start,
            end=end,
            length=end - start,
            detection_method=method_str,
            confidence_score=0.0  # Will be calculated later by ConfidenceScorer
        )
        
        # Combine evidence from both candidates
        # Take the maximum values for coverage-related metrics
        merged.mean_coverage = max(
            c1.mean_coverage or 0, c2.mean_coverage or 0
        ) or None
        
        merged.fold_enrichment = max(
            c1.fold_enrichment or 0, c2.fold_enrichment or 0
        ) or None
        
        merged.coverage_uniformity = max(
            c1.coverage_uniformity or 0, c2.coverage_uniformity or 0
        ) or None
        
        # Sum the support counts
        merged.junction_support = (
            (c1.junction_support or 0) + (c2.junction_support or 0)
        ) or None
        
        merged.split_support = (
            (c1.split_support or 0) + (c2.split_support or 0)
        ) or None
        
        return merged