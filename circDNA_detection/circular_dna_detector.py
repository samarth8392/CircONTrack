#!/usr/bin/env python3
"""
Unified CircDNA Detection Pipeline with Advanced Confidence Scoring
Integrates sophisticated confidence scoring system with optimized detection methods
"""

import sys
import time
import logging
import argparse
from typing import List, Optional, Dict, Any
from tqdm import tqdm
import pysam
import numpy as np
from dataclasses import dataclass
from .coverage_analyzer import CoverageAnalyzer
from .utils import (
    CircularCandidate,
    merge_overlapping_candidates,
    filter_candidates_by_confidence,
    get_top_candidates,
    format_candidate_summary
)


class AdvancedConfidenceScorer:
    """
    Advanced confidence scoring system for CircDNA candidates
    Provides sophisticated scoring with multiple evidence types
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
            'max_split_support': 10.0,
            'min_uniformity': 0.3
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
            # Handle potential negative uniformity scores
            uniformity_score = max(0, min(candidate.coverage_uniformity, 1.0))
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
            # Create temporary candidate for single method scoring
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
            # Two methods: weighted average + bonus
            base_score = np.mean(method_scores)
            bonus = 0.2 * min(method_scores)  # Bonus for multi-method support
            return min(base_score + bonus, 1.0)
        else:
            # Three or more methods: weighted average + larger bonus
            base_score = np.mean(method_scores)
            bonus = 0.3 * min(method_scores)
            return min(base_score + bonus, 1.0)
    
    def get_score_breakdown(self, candidate: CircularCandidate) -> Dict[str, float]:
        """
        Get detailed breakdown of confidence score components
        Useful for debugging and detailed analysis
        """
        breakdown = {}
        
        if candidate.detection_method == 'coverage':
            if candidate.fold_enrichment:
                breakdown['fold_enrichment'] = min(
                    candidate.fold_enrichment / self.thresholds['max_fold_enrichment'], 1.0
                ) * self.weights['coverage']['fold_enrichment']
            
            if candidate.coverage_uniformity:
                breakdown['coverage_uniformity'] = max(0, min(candidate.coverage_uniformity, 1.0)) * \
                    self.weights['coverage']['coverage_uniformity']
            
            if candidate.mean_coverage:
                breakdown['mean_coverage'] = min(
                    candidate.mean_coverage / self.thresholds['max_coverage'], 1.0
                ) * self.weights['coverage']['mean_coverage']
        
        elif candidate.detection_method == 'junction':
            breakdown['base_confidence'] = self.weights['junction']['base_confidence']
            if candidate.junction_support:
                breakdown['junction_support'] = min(
                    candidate.junction_support / self.thresholds['max_junction_support'], 1.0
                ) * self.weights['junction']['junction_support']
        
        elif candidate.detection_method == 'split_read':
            breakdown['base_confidence'] = self.weights['split_read']['base_confidence']
            if candidate.split_support:
                breakdown['split_support'] = min(
                    candidate.split_support / self.thresholds['max_split_support'], 1.0
                ) * self.weights['split_read']['split_support']
        
        return breakdown


class MultiMethodIntegrator:
    """
    Integrates candidates from multiple detection methods with advanced merging
    """
    
    def __init__(self, max_overlap_distance=1000):
        self.max_overlap_distance = max_overlap_distance
        self.confidence_scorer = AdvancedConfidenceScorer()
    
    def integrate_candidates(self, coverage_candidates: List[CircularCandidate],
                           junction_candidates: List[CircularCandidate],
                           split_candidates: List[CircularCandidate]) -> List[CircularCandidate]:
        """
        Integrate candidates from multiple detection methods using utility functions
        """
        all_candidates = coverage_candidates + junction_candidates + split_candidates
        
        if not all_candidates:
            return []
        
        # Sort by chromosome and position
        all_candidates.sort(key=lambda x: (x.chromosome, x.start))
        
        # Use the utility function for merging overlapping candidates
        merged_candidates = merge_overlapping_candidates(all_candidates, self.max_overlap_distance)
        
        # Recalculate confidence scores for merged candidates
        for candidate in merged_candidates:
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        return merged_candidates


class UnifiedCircularDNADetector:
    """
    Unified CircDNA detector with advanced confidence scoring and multi-method integration
    """
    
    def __init__(self, min_fold_enrichment=1.5, min_coverage=5, min_length=200, 
                 max_length=100000, min_confidence=0.3, max_candidates=None, verbose=True):
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.min_confidence = min_confidence
        self.max_candidates = max_candidates  # New parameter for top-N filtering
        self.verbose = verbose
        
        # Initialize advanced components
        self.confidence_scorer = AdvancedConfidenceScorer()
        self.integrator = MultiMethodIntegrator()
        
        # Initialize coverage analyzer
        self.coverage_analyzer = CoverageAnalyzer(
            min_fold_enrichment=min_fold_enrichment,
            min_coverage=min_coverage
        )
        
        # Setup logging
        self._setup_logging()
    
    def _setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO if self.verbose else logging.WARNING,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def detect_circular_dna(self, bam_file: str, reference_file: str, 
                          output_file: Optional[str] = None, 
                          chromosome: Optional[str] = None) -> List[CircularCandidate]:
        """Main detection pipeline with advanced confidence scoring"""
        start_time = time.time()
        self.logger.info("Starting unified CircDNA detection pipeline")
        
        # Validate inputs
        self._validate_inputs(bam_file, reference_file)
        
        # Process chromosomes
        all_candidates = self._process_chromosomes(bam_file, chromosome)
        
        # Integrate and filter results using utility functions
        final_candidates = self._finalize_candidates(all_candidates)
        
        # Write output
        if output_file:
            self._write_output(final_candidates, output_file)
        
        # Print summary
        self._print_summary(final_candidates, time.time() - start_time)
        
        return final_candidates
    
    def _validate_inputs(self, bam_file: str, reference_file: str):
        """Validate input files"""
        self.logger.info("Validating input files")
        
        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                if not bam.has_index():
                    self.logger.warning("BAM file is not indexed - analysis will be slower")
        except Exception as e:
            raise ValueError(f"Invalid BAM file: {str(e)}")
        
        try:
            with pysam.FastaFile(reference_file) as ref:
                self.logger.info(f"Reference validated: {len(ref.references)} contigs")
        except Exception as e:
            raise ValueError(f"Invalid reference file: {str(e)}")
    
    def _process_chromosomes(self, bam_file: str, 
                           target_chromosome: Optional[str] = None) -> List[CircularCandidate]:
        """Process all chromosomes with integrated analysis"""
        self.logger.info("Processing chromosomes")
        
        all_candidates = []
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            chromosomes = [target_chromosome] if target_chromosome else bam.references
            
            for chr_name in tqdm(chromosomes, desc="Chromosomes", disable=not self.verbose):
                try:
                    chr_length = bam.get_reference_length(chr_name)
                    self.logger.info(f"Processing {chr_name} (length: {chr_length:,} bp)")
                    
                    # Run all detection methods for this chromosome
                    coverage_candidates = self._detect_coverage_patterns(bam, chr_name, chr_length)
                    junction_candidates = self._detect_junctions(bam, chr_name)
                    split_candidates = self._detect_split_reads(bam, chr_name)
                    
                    # Integrate results for this chromosome
                    chr_candidates = self.integrator.integrate_candidates(
                        coverage_candidates, junction_candidates, split_candidates
                    )
                    
                    all_candidates.extend(chr_candidates)
                    
                except Exception as e:
                    self.logger.error(f"Error processing {chr_name}: {str(e)}")
                    continue
        
        return all_candidates
    
    def _detect_coverage_patterns(self, bam: pysam.AlignmentFile, 
                                 chromosome: str, chr_length: int) -> List[CircularCandidate]:
        """Detect coverage-based candidates with confidence scoring"""
        
        # Use the optimized coverage analyzer
        coverage_candidates = self.coverage_analyzer.analyze_chromosome_coverage(
            bam, chromosome, chr_length
        )
        
        # Convert to CircularCandidate objects with confidence scores
        candidates = []
        for cov_candidate in coverage_candidates:
            candidate = CircularCandidate(
                chromosome=cov_candidate['chromosome'],
                start=cov_candidate['start'],
                end=cov_candidate['end'],
                length=cov_candidate['length'],
                detection_method='coverage',
                mean_coverage=cov_candidate['mean_coverage'],
                fold_enrichment=cov_candidate['fold_enrichment'],
                coverage_uniformity=cov_candidate['coverage_uniformity']
            )
            
            # Calculate confidence score using advanced scorer
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
            candidates.append(candidate)
        
        return candidates

    def _detect_junctions(self, bam: pysam.AlignmentFile, 
                         chromosome: str) -> List[CircularCandidate]:
        """Detect junction-based candidates with confidence scoring"""
        candidates = []
        junction_positions = {}
        
        # Collect junction evidence
        for read in bam.fetch(chromosome):
            if self._is_junction_read(read):
                pos = read.reference_start
                junction_positions[pos] = junction_positions.get(pos, 0) + 1
        
        # Create candidates from junction clusters
        for pos, support in junction_positions.items():
            if support >= 3:  # Minimum support threshold
                candidate = CircularCandidate(
                    chromosome=chromosome,
                    start=pos - 500,  # Estimate boundaries
                    end=pos + 500,
                    length=1000,
                    detection_method='junction',
                    junction_support=support
                )
                
                # Calculate confidence score using advanced scorer
                candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
                candidates.append(candidate)
        
        return candidates
    
    def _detect_split_reads(self, bam: pysam.AlignmentFile, 
                           chromosome: str) -> List[CircularCandidate]:
        """Detect split-read based candidates with confidence scoring"""
        candidates = []
        split_positions = {}
        
        # Collect split read evidence
        for read in bam.fetch(chromosome):
            if self._is_split_read(read):
                pos = read.reference_start
                split_positions[pos] = split_positions.get(pos, 0) + 1
        
        # Create candidates from split clusters
        for pos, support in split_positions.items():
            if support >= 3:  # Minimum support threshold
                candidate = CircularCandidate(
                    chromosome=chromosome,
                    start=pos - 500,
                    end=pos + 500,
                    length=1000,
                    detection_method='split_read',
                    split_support=support
                )
                
                # Calculate confidence score using advanced scorer
                candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
                candidates.append(candidate)
        
        return candidates
    
    def _finalize_candidates(self, candidates: List[CircularCandidate]) -> List[CircularCandidate]:
        """Final filtering and sorting of candidates using utility functions"""
        # Filter by length
        length_filtered = [c for c in candidates 
                          if self.min_length <= c.length <= self.max_length]
        
        # Use utility function for confidence filtering
        confidence_filtered = filter_candidates_by_confidence(length_filtered, self.min_confidence)
        
        # Sort by confidence score (descending)
        sorted_candidates = sorted(confidence_filtered, key=lambda x: x.confidence_score, reverse=True)
        
        # Apply top-N filtering if specified
        if self.max_candidates is not None:
            return get_top_candidates(sorted_candidates, self.max_candidates)
        
        return sorted_candidates
    
    # Helper methods for read classification
    def _is_junction_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read shows junction signature"""
        return (not read.is_unmapped and 
                read.mapping_quality >= 20 and
                read.cigarstring and 'S' in read.cigarstring)
    
    def _is_split_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read is split alignment"""
        return (read.has_tag('SA') or 
                (read.cigarstring and read.cigarstring.count('S') >= 2))
    
    # Output methods
    def _write_output(self, candidates: List[CircularCandidate], output_file: str):
        """Write results to BED file with confidence scores"""
        self.logger.info(f"Writing {len(candidates)} candidates to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("# Unified CircDNA Detection Results with Advanced Confidence Scoring\n")
            f.write("# chr\tstart\tend\tname\tconfidence\tstrand\tmethod\tlength\tevidence\n")
            
            for i, candidate in enumerate(candidates):
                # Format evidence string
                evidence_parts = []
                if candidate.fold_enrichment:
                    evidence_parts.append(f"fold={candidate.fold_enrichment:.2f}")
                if candidate.mean_coverage:
                    evidence_parts.append(f"cov={candidate.mean_coverage:.1f}")
                if candidate.junction_support:
                    evidence_parts.append(f"junc={candidate.junction_support}")
                if candidate.split_support:
                    evidence_parts.append(f"split={candidate.split_support}")
                
                evidence_str = ";".join(evidence_parts) if evidence_parts else "none"
                
                f.write(f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
                       f"circDNA_{i+1}\t{candidate.confidence_score:.3f}\t.\t"
                       f"{candidate.detection_method}\t{candidate.length}\t{evidence_str}\n")
    
    def _print_summary(self, candidates: List[CircularCandidate], elapsed_time: float):
        """Print comprehensive analysis summary using utility functions"""
        print("\n" + "="*70)
        print("UNIFIED CIRCULAR DNA DETECTION SUMMARY")
        print("="*70)
        print(f"Total high-confidence candidates: {len(candidates)}")
        print(f"Analysis time: {elapsed_time:.2f} seconds")
        print(f"Confidence threshold: {self.min_confidence:.3f}")
        
        if candidates:
            # Method breakdown
            method_counts = {}
            for candidate in candidates:
                method = candidate.detection_method
                method_counts[method] = method_counts.get(method, 0) + 1
            
            print("\nCandidates by detection method:")
            for method, count in sorted(method_counts.items()):
                print(f"  {method}: {count}")
            
            # Confidence distribution
            confidences = [c.confidence_score for c in candidates]
            print(f"\nConfidence score distribution:")
            print(f"  Mean: {np.mean(confidences):.3f}")
            print(f"  Median: {np.median(confidences):.3f}")
            print(f"  Range: {min(confidences):.3f} - {max(confidences):.3f}")
            
            # Top candidates with detailed information using utility function
            top_candidates = get_top_candidates(candidates, 5)
            print(f"\nTop 5 candidates with detailed scoring:")
            for i, candidate in enumerate(top_candidates):
                print(f"\n  {i+1}. {format_candidate_summary(candidate)}")
                print(f"     Confidence: {candidate.confidence_score:.3f}")
                print(f"     Length: {candidate.length:,} bp")
                
                # Show confidence score breakdown
                breakdown = self.confidence_scorer.get_score_breakdown(candidate)
                if breakdown:
                    print(f"     Score breakdown: {breakdown}")
        
        print("="*70)
    
    def get_detailed_report(self, candidates: List[CircularCandidate]) -> str:
        """Generate detailed report for candidates using utility functions"""
        report = []
        report.append("DETAILED CIRCULAR DNA DETECTION REPORT")
        report.append("=" * 50)
        
        for i, candidate in enumerate(candidates):
            report.append(f"\nCandidate {i+1}:")
            report.append(f"  Summary: {format_candidate_summary(candidate)}")
            report.append(f"  Confidence score: {candidate.confidence_score:.3f}")
            
            # Detailed evidence
            if candidate.fold_enrichment:
                report.append(f"  Fold enrichment: {candidate.fold_enrichment:.2f}")
            if candidate.mean_coverage:
                report.append(f"  Mean coverage: {candidate.mean_coverage:.1f}")
            if candidate.coverage_uniformity:
                report.append(f"  Coverage uniformity: {candidate.coverage_uniformity:.3f}")
            if candidate.junction_support:
                report.append(f"  Junction support: {candidate.junction_support}")
            if candidate.split_support:
                report.append(f"  Split read support: {candidate.split_support}")
            
            # Score breakdown
            breakdown = self.confidence_scorer.get_score_breakdown(candidate)
            if breakdown:
                report.append(f"  Score breakdown:")
                for component, score in breakdown.items():
                    report.append(f"    {component}: {score:.3f}")
        
        return "\n".join(report)


def main():
    """Main function with comprehensive argument parsing"""
    parser = argparse.ArgumentParser(
        description='Unified CircDNA Detection with Advanced Confidence Scoring',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.bam reference.fa -o results.bed
  %(prog)s input.bam reference.fa -c chr1 --min-confidence 0.5 --max-candidates 20
  %(prog)s input.bam reference.fa --detailed-report results_detailed.txt
        """
    )
    
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('reference_file', help='Reference FASTA file')
    parser.add_argument('-o', '--output', default='circular_dna_unified.bed', 
                       help='Output BED file (default: circular_dna_unified.bed)')
    parser.add_argument('-c', '--chromosome', help='Analyze specific chromosome only')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--detailed-report', help='Generate detailed report file')
    
    # Detection parameters
    parser.add_argument('--min-fold-enrichment', type=float, default=1.5,
                       help='Minimum fold enrichment (default: 1.5)')
    parser.add_argument('--min-coverage', type=int, default=5,
                       help='Minimum coverage depth (default: 5)')
    parser.add_argument('--min-length', type=int, default=200,
                       help='Minimum circular DNA length (default: 200)')
    parser.add_argument('--max-length', type=int, default=100000,
                       help='Maximum circular DNA length (default: 100000)')
    parser.add_argument('--min-confidence', type=float, default=0.3,
                       help='Minimum confidence score (default: 0.3)')
    parser.add_argument('--max-candidates', type=int, 
                       help='Maximum number of candidates to report (top N by confidence)')

    args = parser.parse_args()
    
    # Create unified detector
    detector = UnifiedCircularDNADetector(
        min_fold_enrichment=args.min_fold_enrichment,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        max_length=args.max_length,
        min_confidence=args.min_confidence,
        max_candidates=args.max_candidates,
        verbose=not args.quiet
    )
    
    # Run detection
    results = detector.detect_circular_dna(
        args.bam_file,
        args.reference_file,
        args.output,
        args.chromosome
    )
    
    # Generate detailed report if requested
    if args.detailed_report:
        detailed_report = detector.get_detailed_report(results)
        with open(args.detailed_report, 'w') as f:
            f.write(detailed_report)
        print(f"Detailed report written to {args.detailed_report}")
    
    print(f"\nAnalysis complete! Found {len(results)} high-confidence candidates.")
    if args.output:
        print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()