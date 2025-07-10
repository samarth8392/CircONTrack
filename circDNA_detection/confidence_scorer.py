#!/usr/bin/env python3
"""
CircDNA Pipeline with Proper Confidence Scoring System
Fixes the missing confidence scoring issue identified in the analysis
"""

import numpy as np
from typing import List, Dict, Tuple
from dataclasses import dataclass
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
        
        # Combine scores with diminishing returns
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


# Updated detection modules with confidence scoring

class UpdatedCoverageAnalyzer:
    """Coverage analyzer with confidence scoring"""
    
    def __init__(self, window_sizes=[100, 500, 1000], min_fold_enrichment=1.5, 
                 min_coverage=5, uniformity_threshold=0.4):
        self.window_sizes = window_sizes
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.uniformity_threshold = uniformity_threshold
        self.confidence_scorer = ConfidenceScorer()
    
    def _create_candidate_with_confidence(self, chromosome, start, end, length, 
                                        mean_coverage, fold_enrichment, coverage_uniformity):
        """Create candidate with proper confidence scoring"""
        candidate = CircularCandidate(
            chromosome=chromosome,
            start=start,
            end=end,
            length=length,
            mean_coverage=mean_coverage,
            fold_enrichment=fold_enrichment,
            coverage_uniformity=coverage_uniformity,
            confidence_score=0.0,  # Will be calculated below
            detection_method='coverage'
        )
        
        # Calculate confidence score
        candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        return candidate


class UpdatedJunctionDetector:
    """Junction detector with confidence scoring"""
    
    def __init__(self, min_support=3, max_junction_distance=1000):
        self.min_support = min_support
        self.max_junction_distance = max_junction_distance
        self.confidence_scorer = ConfidenceScorer()
    
    def _create_candidate_with_confidence(self, chromosome, start, end, length, junction_support):
        """Create candidate with proper confidence scoring"""
        candidate = CircularCandidate(
            chromosome=chromosome,
            start=start,
            end=end,
            length=length,
            junction_support=junction_support,
            confidence_score=0.0,  # Will be calculated below
            detection_method='junction'
        )
        
        # Calculate confidence score
        candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        return candidate


class UpdatedSplitReadAnalyzer:
    """Split read analyzer with confidence scoring"""
    
    def __init__(self, min_split_length=50, min_support=3, max_distance=1000):
        self.min_split_length = min_split_length
        self.min_support = min_support
        self.max_distance = max_distance
        self.confidence_scorer = ConfidenceScorer()
    
    def _create_candidate_with_confidence(self, chromosome, start, end, length, split_support):
        """Create candidate with proper confidence scoring"""
        candidate = CircularCandidate(
            chromosome=chromosome,
            start=start,
            end=end,
            length=length,
            split_support=split_support,
            confidence_score=0.0,  # Will be calculated below
            detection_method='split_read'
        )
        
        # Calculate confidence score
        candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        return candidate


class MultiMethodIntegrator:
    """
    Integrates candidates from multiple detection methods
    """
    
    def __init__(self, max_overlap_distance=1000):
        self.max_overlap_distance = max_overlap_distance
        self.confidence_scorer = ConfidenceScorer()
    
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
        
        # Recalculate confidence scores for merged candidates
        for candidate in merged_candidates:
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        return merged_candidates
    
    def _merge_overlapping_candidates(self, candidates: List[CircularCandidate]) -> List[CircularCandidate]:
        """Merge overlapping candidates from different methods"""
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
            confidence_score=0.0  # Will be calculated later
        )
        
        # Combine evidence from both candidates
        merged.mean_coverage = max(
            c1.mean_coverage or 0, c2.mean_coverage or 0
        ) or None
        
        merged.fold_enrichment = max(
            c1.fold_enrichment or 0, c2.fold_enrichment or 0
        ) or None
        
        merged.coverage_uniformity = max(
            c1.coverage_uniformity or 0, c2.coverage_uniformity or 0
        ) or None
        
        merged.junction_support = (
            (c1.junction_support or 0) + (c2.junction_support or 0)
        ) or None
        
        merged.split_support = (
            (c1.split_support or 0) + (c2.split_support or 0)
        ) or None
        
        return merged


class UpdatedCircularDNADetector:
    """
    Main detector with proper confidence scoring and multi-method integration
    """
    
    def __init__(self, min_fold_enrichment=1.5, min_coverage=5, min_length=200, 
                 max_length=100000, min_confidence=0.3, verbose=False):
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.min_confidence = min_confidence
        self.verbose = verbose
        
        # Initialize analyzers
        self.coverage_analyzer = UpdatedCoverageAnalyzer(
            min_fold_enrichment=min_fold_enrichment,
            min_coverage=min_coverage
        )
        self.junction_detector = UpdatedJunctionDetector()
        self.split_analyzer = UpdatedSplitReadAnalyzer()
        self.integrator = MultiMethodIntegrator()
    
    def detect_circular_dna(self, bam_file, reference_file, output_file=None, 
                          chromosome=None):
        """
        Main detection pipeline with confidence scoring
        """
        if self.verbose:
            print("Starting CircDNA detection with confidence scoring...")
        
        # Run individual detection methods
        coverage_candidates = self.coverage_analyzer.detect_coverage_patterns(
            bam_file, chromosome
        )
        
        junction_candidates = self.junction_detector.detect_junctions(
            bam_file, chromosome
        )
        
        split_candidates = self.split_analyzer.analyze_split_reads(
            bam_file, chromosome
        )
        
        # Integrate results
        integrated_candidates = self.integrator.integrate_candidates(
            coverage_candidates, junction_candidates, split_candidates
        )
        
        # Filter by confidence threshold
        high_confidence_candidates = [
            c for c in integrated_candidates 
            if c.confidence_score >= self.min_confidence
        ]
        
        # Sort by confidence score (descending)
        final_candidates = sorted(
            high_confidence_candidates, 
            key=lambda x: x.confidence_score, 
            reverse=True
        )
        
        if output_file:
            self._write_output_with_confidence(final_candidates, output_file)
        
        if self.verbose:
            self._print_confidence_summary(
                coverage_candidates, junction_candidates, split_candidates, 
                integrated_candidates, final_candidates
            )
        
        return final_candidates
    
    def _write_output_with_confidence(self, candidates, output_file):
        """Write results with proper confidence scores"""
        with open(output_file, 'w') as f:
            f.write("# CircDNA Detection Results with Confidence Scores\n")
            f.write("# chr\tstart\tend\tname\tconfidence\tstrand\tmethod\tlength\n")
            
            for i, candidate in enumerate(candidates):
                f.write(f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
                       f"circDNA_{i+1}\t{candidate.confidence_score:.3f}\t.\t"
                       f"{candidate.detection_method}\t{candidate.length}\n")
    
    def _print_confidence_summary(self, coverage_cands, junction_cands, split_cands, 
                                integrated_cands, final_cands):
        """Print comprehensive confidence scoring summary"""
        print("\n" + "="*60)
        print("CONFIDENCE SCORING SUMMARY")
        print("="*60)
        
        print(f"Raw candidates found:")
        print(f"  Coverage-based: {len(coverage_cands)}")
        print(f"  Junction-based: {len(junction_cands)}")
        print(f"  Split-read based: {len(split_cands)}")
        print(f"  Total raw: {len(coverage_cands) + len(junction_cands) + len(split_cands)}")
        
        print(f"\nAfter multi-method integration: {len(integrated_cands)}")
        print(f"After confidence filtering (≥{self.min_confidence:.2f}): {len(final_cands)}")
        
        if final_cands:
            print(f"\nTop 5 candidates by confidence:")
            for i, candidate in enumerate(final_cands[:5]):
                print(f"  {i+1}. {candidate.chromosome}:{candidate.start}-{candidate.end}")
                print(f"     Confidence: {candidate.confidence_score:.3f}")
                print(f"     Method: {candidate.detection_method}")
                print(f"     Length: {candidate.length:,} bp")
                
                # Show contributing evidence
                evidence = []
                if candidate.fold_enrichment:
                    evidence.append(f"fold={candidate.fold_enrichment:.1f}")
                if candidate.junction_support:
                    evidence.append(f"junctions={candidate.junction_support}")
                if candidate.split_support:
                    evidence.append(f"splits={candidate.split_support}")
                
                if evidence:
                    print(f"     Evidence: {', '.join(evidence)}")
                print()
        
        print("="*60)


# Example usage and testing
def test_confidence_scoring():
    """Test the confidence scoring system"""
    scorer = ConfidenceScorer()
    
    # Test coverage-based candidate
    coverage_candidate = CircularCandidate(
        chromosome='chr1',
        start=1000,
        end=2000,
        length=1000,
        mean_coverage=15.0,
        fold_enrichment=3.2,
        coverage_uniformity=0.8,
        detection_method='coverage'
    )
    
    coverage_score = scorer.calculate_confidence(coverage_candidate)
    print(f"Coverage candidate confidence: {coverage_score:.3f}")
    
    # Test junction-based candidate
    junction_candidate = CircularCandidate(
        chromosome='chr1',
        start=1000,
        end=2000,
        length=1000,
        junction_support=8,
        detection_method='junction'
    )
    
    junction_score = scorer.calculate_confidence(junction_candidate)
    print(f"Junction candidate confidence: {junction_score:.3f}")
    
    # Test multi-method candidate
    multi_candidate = CircularCandidate(
        chromosome='chr1',
        start=1000,
        end=2000,
        length=1000,
        mean_coverage=12.0,
        fold_enrichment=2.5,
        coverage_uniformity=0.7,
        junction_support=5,
        split_support=3,
        detection_method='coverage+junction+split_read'
    )
    
    multi_score = scorer.calculate_confidence(multi_candidate)
    print(f"Multi-method candidate confidence: {multi_score:.3f}")


if __name__ == "__main__":
    test_confidence_scoring()