#!/usr/bin/env python3
"""
FIXED CircONTrack Circular DNA Detector
Integrates the sophisticated analysis modules that were missing from the main detector
"""
import sys
import time
import logging
from tqdm import tqdm
import pysam
import numpy as np
from scipy.stats import median_abs_deviation
from scipy.signal import find_peaks
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Union

@dataclass
class CircularCandidate:
    """Data class for circular DNA candidates"""
    chromosome: str
    start: int
    end: int
    length: int
    detection_method: str
    confidence_score: float = 0.0
    
    # Coverage-based evidence
    mean_coverage: Optional[float] = None
    fold_enrichment: Optional[float] = None
    coverage_uniformity: Optional[float] = None
    
    # Junction-based evidence
    junction_support: Optional[int] = None
    
    # Split-read evidence
    split_support: Optional[int] = None

class ConfidenceScorer:
    """Confidence scoring system for CircDNA candidates"""
    
    def __init__(self):
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
        
        self.thresholds = {
            'max_fold_enrichment': 5.0,
            'max_coverage': 20.0,
            'max_junction_support': 10.0,
            'max_split_support': 10.0
        }
    
    def calculate_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence score based on available evidence"""
        if candidate.detection_method == 'coverage':
            return self._calculate_coverage_confidence(candidate)
        elif candidate.detection_method == 'junction':
            return self._calculate_junction_confidence(candidate)
        elif candidate.detection_method == 'split_read':
            return self._calculate_split_read_confidence(candidate)
        elif '+' in candidate.detection_method:
            return self._calculate_multi_method_confidence(candidate)
        else:
            return 0.1
    
    def _calculate_coverage_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for coverage-based detection"""
        score = 0.0
        weights = self.weights['coverage']
        
        # Fold enrichment component (CRITICAL)
        if candidate.fold_enrichment is not None:
            fold_score = min(candidate.fold_enrichment / self.thresholds['max_fold_enrichment'], 1.0)
            score += fold_score * weights['fold_enrichment']
        
        # Coverage uniformity component
        if candidate.coverage_uniformity is not None:
            uniformity_score = max(0, candidate.coverage_uniformity)
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
        
        if candidate.junction_support is not None:
            support_score = min(candidate.junction_support / self.thresholds['max_junction_support'], 1.0)
            score += support_score * weights['junction_support']
        
        score += weights['base_confidence']
        return min(score, 1.0)
    
    def _calculate_split_read_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for split-read detection"""
        score = 0.0
        weights = self.weights['split_read']
        
        if candidate.split_support is not None:
            support_score = min(candidate.split_support / self.thresholds['max_split_support'], 1.0)
            score += support_score * weights['split_support']
        
        score += weights['base_confidence']
        return min(score, 1.0)
    
    def _calculate_multi_method_confidence(self, candidate: CircularCandidate) -> float:
        """Calculate confidence for multi-method detection"""
        methods = candidate.detection_method.split('+')
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
        
        if len(method_scores) == 1:
            return method_scores[0]
        elif len(method_scores) == 2:
            base_score = np.mean(method_scores)
            bonus = 0.2 * min(method_scores)
            return min(base_score + bonus, 1.0)
        else:
            base_score = np.mean(method_scores)
            bonus = 0.3 * min(method_scores)
            return min(base_score + bonus, 1.0)

class ProperCoverageAnalyzer:
    """FIXED Coverage analyzer with proper statistical analysis"""
    
    def __init__(self, window_sizes=[100, 500, 1000], min_fold_enrichment=1.5, 
                 min_coverage=5, uniformity_threshold=0.4):
        self.window_sizes = window_sizes
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.uniformity_threshold = uniformity_threshold
        self.confidence_scorer = ConfidenceScorer()
    
    def analyze_coverage(self, bam, chromosome, verbose=False):
        """PROPER coverage analysis with statistical methods"""
        if verbose:
            print(f"  → Analyzing coverage for {chromosome} using proper statistics")
        
        all_candidates = []
        
        # Multi-scale analysis
        for window_size in self.window_sizes:
            candidates = self._analyze_chromosome_coverage(
                bam, chromosome, window_size, verbose
            )
            all_candidates.extend(candidates)
        
        # Consolidate overlapping candidates
        consolidated = self._consolidate_candidates(all_candidates)
        
        if verbose:
            print(f"  → Coverage analysis complete: {len(consolidated)} high-confidence candidates")
        
        return consolidated
    
    def _analyze_chromosome_coverage(self, bam, chromosome, window_size, verbose=False):
        """Analyze coverage patterns for single chromosome at given window size"""
        chrom_length = bam.get_reference_length(chromosome)
        if chrom_length < window_size * 10:
            return []
        
        # Calculate coverage in windows
        coverage_data = self._calculate_windowed_coverage(
            bam, chromosome, window_size, chrom_length, verbose
        )
        
        if len(coverage_data) < 10:
            return []
        
        # PROPER STATISTICAL ANALYSIS
        coverage_values = np.array([w['coverage'] for w in coverage_data])
        positions = np.array([w['start'] for w in coverage_data])
        
        # Use median and MAD for robust threshold calculation
        median_cov = np.median(coverage_values)
        mad_cov = median_abs_deviation(coverage_values)
        
        if median_cov < 1 or mad_cov == 0:
            return []
        
        # CRITICAL: Proper fold enrichment threshold
        fold_threshold = median_cov * self.min_fold_enrichment
        robust_threshold = median_cov + 3 * mad_cov
        threshold = max(fold_threshold, robust_threshold, self.min_coverage)
        
        if verbose:
            print(f"    Window size {window_size}: median={median_cov:.1f}, "
                  f"MAD={mad_cov:.1f}, threshold={threshold:.1f}")
        
        # Find peaks with proper background subtraction
        candidates = self._find_coverage_peaks(
            coverage_values, positions, threshold, median_cov, window_size, chromosome
        )
        
        return candidates
    
    def _calculate_windowed_coverage(self, bam, chromosome, window_size, chrom_length, verbose):
        """Calculate coverage in sliding windows with quality filtering"""
        coverage_data = []
        
        progress_iter = range(0, chrom_length, window_size)
        if verbose:
            progress_iter = tqdm(progress_iter, desc=f"    Windows ({window_size}bp)", 
                               leave=False, file=sys.stdout)
        
        for start in progress_iter:
            end = min(start + window_size, chrom_length)
            
            # Count high-quality reads only
            count = 0
            for read in bam.fetch(chromosome, start, end):
                if (read.mapping_quality >= 20 and 
                    not read.is_secondary and 
                    not read.is_duplicate and
                    not read.is_unmapped):
                    count += 1
            
            # Normalize by window size
            normalized_coverage = count / (end - start) * window_size
            
            coverage_data.append({
                'start': start,
                'end': end,
                'coverage': normalized_coverage,
                'raw_count': count
            })
        
        return coverage_data
    
    def _find_coverage_peaks(self, coverage_values, positions, threshold, baseline, window_size, chromosome):
        """Find coverage peaks using proper peak detection algorithms"""
        candidates = []
        
        # LOCAL BACKGROUND SUBTRACTION
        local_bg = self._calculate_local_background(coverage_values, window_size=10)
        adjusted_coverage = coverage_values - local_bg
        
        # PROPER PEAK DETECTION
        peak_indices, peak_properties = find_peaks(
            adjusted_coverage,
            height=threshold - baseline,
            distance=max(5, 500 // window_size),
            width=max(2, 200 // window_size),
            prominence=baseline * 0.5  # Require significant prominence
        )
        
        for peak_idx in peak_indices:
            # Define peak boundaries
            start_idx, end_idx = self._find_peak_boundaries(
                adjusted_coverage, peak_idx, baseline * 0.5
            )
            
            # Convert to genomic coordinates
            genomic_start = positions[start_idx]
            genomic_end = positions[min(end_idx, len(positions)-1)] + window_size
            region_length = genomic_end - genomic_start
            
            # Filter by length requirements
            if region_length < 200 or region_length > 100000:
                continue
            
            # Calculate region statistics
            region_coverage = coverage_values[start_idx:end_idx+1]
            region_median = np.median(region_coverage)
            region_mad = median_abs_deviation(region_coverage)
            
            # Calculate flanking background
            flanking_bg = self._calculate_flanking_background(
                coverage_values, start_idx, end_idx, baseline
            )
            
            # CRITICAL: Calculate actual fold enrichment
            fold_enrichment = region_median / max(flanking_bg, 1)
            
            # Only keep candidates with sufficient fold enrichment
            if fold_enrichment < self.min_fold_enrichment:
                continue
            
            # Calculate uniformity score
            uniformity = 1 - (region_mad / max(region_median, 1))
            
            # Create candidate with proper confidence scoring
            candidate = CircularCandidate(
                chromosome=chromosome,
                start=genomic_start,
                end=genomic_end,
                length=region_length,
                mean_coverage=region_median,
                fold_enrichment=fold_enrichment,
                coverage_uniformity=uniformity,
                detection_method='coverage'
            )
            
            # Calculate confidence score
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
            
            candidates.append(candidate)
        
        return candidates
    
    def _calculate_local_background(self, coverage_values, window_size=10):
        """Calculate local background using sliding median"""
        local_bg = np.zeros_like(coverage_values)
        
        for i in range(len(coverage_values)):
            start = max(0, i - window_size)
            end = min(len(coverage_values), i + window_size + 1)
            local_bg[i] = np.median(coverage_values[start:end])
        
        return local_bg
    
    def _find_peak_boundaries(self, coverage, peak_idx, threshold):
        """Find peak boundaries using threshold crossing"""
        # Left boundary
        left_boundary = peak_idx
        for i in range(peak_idx - 1, -1, -1):
            if coverage[i] <= threshold:
                left_boundary = i + 1
                break
            left_boundary = i
        
        # Right boundary
        right_boundary = peak_idx
        for i in range(peak_idx + 1, len(coverage)):
            if coverage[i] <= threshold:
                right_boundary = i - 1
                break
            right_boundary = i
        
        return left_boundary, right_boundary
    
    def _calculate_flanking_background(self, coverage, start_idx, end_idx, baseline):
        """Calculate background coverage from flanking regions"""
        region_size = end_idx - start_idx
        flank_size = min(region_size, 20)
        
        # Left flank
        left_start = max(0, start_idx - flank_size)
        left_flank = coverage[left_start:start_idx]
        
        # Right flank
        right_end = min(len(coverage), end_idx + flank_size)
        right_flank = coverage[end_idx:right_end]
        
        # Combine flanks
        flanking_values = np.concatenate([left_flank, right_flank])
        
        return np.median(flanking_values) if len(flanking_values) > 0 else baseline
    
    def _consolidate_candidates(self, candidates):
        """Consolidate overlapping candidates from multi-scale analysis"""
        if not candidates:
            return []
        
        # Sort by position
        candidates.sort(key=lambda x: (x.chromosome, x.start))
        
        consolidated = []
        current = candidates[0]
        
        for candidate in candidates[1:]:
            # Check for overlap
            if (candidate.chromosome == current.chromosome and 
                candidate.start <= current.end + 1000):
                
                # Keep candidate with higher confidence
                if candidate.confidence_score > current.confidence_score:
                    current = candidate
            else:
                consolidated.append(current)
                current = candidate
        
        consolidated.append(current)
        return consolidated

class FixedCircularDNADetector:
    """FIXED CircDNA detector using proper statistical analysis"""
    
    def __init__(self, min_fold_enrichment=1.5, min_coverage=5, min_length=200, 
                 max_length=100000, min_confidence=0.3, verbose=True):
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.min_confidence = min_confidence
        self.verbose = verbose
        
        # Initialize PROPER analyzers
        self.coverage_analyzer = ProperCoverageAnalyzer(
            min_fold_enrichment=min_fold_enrichment,
            min_coverage=min_coverage
        )
        
        self.confidence_scorer = ConfidenceScorer()
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO if verbose else logging.WARNING,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def detect_circular_dna(self, bam_file, reference_file, output_file=None, 
                          chromosome=None):
        """FIXED detection pipeline with proper statistical analysis"""
        start_time = time.time()
        
        if self.verbose:
            print("="*60)
            print("FIXED CircDNA Detection Pipeline")
            print("="*60)
            print(f"Parameters:")
            print(f"  Min fold enrichment: {self.min_fold_enrichment}")
            print(f"  Min coverage: {self.min_coverage}")
            print(f"  Min length: {self.min_length}")
            print(f"  Max length: {self.max_length}")
            print(f"  Min confidence: {self.min_confidence}")
            print()
        
        # Validate inputs
        self._validate_inputs(bam_file, reference_file)
        
        # Process chromosomes
        all_candidates = []
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chr_name in chromosomes:
                if self.verbose:
                    print(f"Processing chromosome {chr_name}...")
                
                # PROPER coverage analysis
                coverage_candidates = self.coverage_analyzer.analyze_coverage(
                    bam, chr_name, self.verbose
                )
                
                # Filter by length and confidence
                filtered_candidates = []
                for candidate in coverage_candidates:
                    if (self.min_length <= candidate.length <= self.max_length and
                        candidate.confidence_score >= self.min_confidence):
                        filtered_candidates.append(candidate)
                
                all_candidates.extend(filtered_candidates)
                
                if self.verbose:
                    print(f"  → {len(coverage_candidates)} raw candidates")
                    print(f"  → {len(filtered_candidates)} high-confidence candidates")
        
        # Sort by confidence score
        final_candidates = sorted(all_candidates, 
                                key=lambda x: x.confidence_score, 
                                reverse=True)
        
        # Write output
        if output_file:
            self._write_output(final_candidates, output_file)
        
        # Print summary
        self._print_summary(final_candidates, start_time)
        
        return final_candidates
    
    def _validate_inputs(self, bam_file, reference_file):
        """Validate input files"""
        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                if not bam.has_index():
                    self.logger.warning("BAM file is not indexed - performance will be slow")
        except Exception as e:
            raise ValueError(f"Invalid BAM file: {str(e)}")
        
        try:
            with pysam.FastaFile(reference_file) as ref:
                pass
        except Exception as e:
            raise ValueError(f"Invalid reference file: {str(e)}")
    
    def _write_output(self, candidates, output_file):
        """Write results with confidence scores"""
        if self.verbose:
            print(f"Writing {len(candidates)} candidates to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("# FIXED CircDNA Detection Results\n")
            f.write("# chr\tstart\tend\tname\tconfidence\tfold_enrichment\tcoverage\tlength\n")
            
            for i, candidate in enumerate(candidates):
                f.write(f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
                       f"circDNA_{i+1}\t{candidate.confidence_score:.3f}\t"
                       f"{candidate.fold_enrichment:.2f}\t{candidate.mean_coverage:.1f}\t"
                       f"{candidate.length}\n")
    
    def _print_summary(self, candidates, start_time):
        """Print comprehensive summary"""
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*60)
        print("FIXED CircDNA DETECTION SUMMARY")
        print("="*60)
        print(f"Total high-confidence candidates: {len(candidates)}")
        print(f"Analysis completed in: {elapsed_time:.2f} seconds")
        
        if candidates:
            # Confidence distribution
            confidences = [c.confidence_score for c in candidates]
            fold_enrichments = [c.fold_enrichment for c in candidates if c.fold_enrichment]
            
            print(f"\nConfidence score distribution:")
            print(f"  Mean: {np.mean(confidences):.3f}")
            print(f"  Median: {np.median(confidences):.3f}")
            print(f"  Min: {min(confidences):.3f}")
            print(f"  Max: {max(confidences):.3f}")
            
            if fold_enrichments:
                print(f"\nFold enrichment distribution:")
                print(f"  Mean: {np.mean(fold_enrichments):.2f}")
                print(f"  Median: {np.median(fold_enrichments):.2f}")
                print(f"  Min: {min(fold_enrichments):.2f}")
                print(f"  Max: {max(fold_enrichments):.2f}")
            
            # Top candidates
            print(f"\nTop 10 candidates by confidence:")
            for i, candidate in enumerate(candidates[:10]):
                print(f"  {i+1}. {candidate.chromosome}:{candidate.start}-{candidate.end}")
                print(f"     Confidence: {candidate.confidence_score:.3f}")
                print(f"     Fold enrichment: {candidate.fold_enrichment:.2f}")
                print(f"     Coverage: {candidate.mean_coverage:.1f}")
                print(f"     Length: {candidate.length:,} bp")
        
        print("="*60)

def main():
    """Main function with proper argument parsing"""
    import argparse
    
    parser = argparse.ArgumentParser(description='FIXED CircDNA Detection Pipeline')
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('reference_file', help='Reference FASTA file')
    parser.add_argument('-o', '--output', help='Output file', default='fixed_circular_dna.bed')
    parser.add_argument('-c', '--chromosome', help='Analyze specific chromosome only')
    parser.add_argument('--min-fold-enrichment', type=float, default=1.5, 
                       help='Minimum fold enrichment (default: 1.5)')
    parser.add_argument('--min-coverage', type=int, default=5,
                       help='Minimum coverage (default: 5)')
    parser.add_argument('--min-length', type=int, default=200,
                       help='Minimum length (default: 200)')
    parser.add_argument('--max-length', type=int, default=100000,
                       help='Maximum length (default: 100000)')
    parser.add_argument('--min-confidence', type=float, default=0.3,
                       help='Minimum confidence score (default: 0.3)')
    parser.add_argument('--quiet', action='store_true', 
                       help='Disable verbose output')

    args = parser.parse_args()
    
    # Create FIXED detector
    detector = FixedCircularDNADetector(
        min_fold_enrichment=args.min_fold_enrichment,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        max_length=args.max_length,
        min_confidence=args.min_confidence,
        verbose=not args.quiet
    )
    
    # Run detection
    results = detector.detect_circular_dna(
        args.bam_file,
        args.reference_file,
        args.output,
        args.chromosome
    )
    
    print(f"\nFixed analysis complete! Results written to {args.output}")
    print(f"Found {len(results)} high-confidence circular DNA candidates")

if __name__ == "__main__":
    main()