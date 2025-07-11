#!/usr/bin/env python3
"""
ONT-Optimized Coverage Pattern Analysis
Multi-scale window analysis with robust statistics
"""

import pysam
import numpy as np
from collections import defaultdict
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation
from tqdm import tqdm
import sys
from .utils import CircularCandidate


class CoverageAnalyzer:
    """Coverage-based circular DNA detection using multi-scale analysis"""
    
    def __init__(self, window_sizes=[500, 1000, 5000], min_fold_enrichment=1.5, 
                 min_coverage=5, uniformity_threshold=0.4, verbose=True):
        self.window_sizes = window_sizes
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.uniformity_threshold = uniformity_threshold
        self.verbose = verbose
    
    def detect_coverage_patterns(self, bam_file, chromosome=None):
        """Multi-scale coverage pattern detection"""
        all_candidates = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            # Create chromosome progress bar
            chr_progress = tqdm(chromosomes, 
                              desc="Coverage analysis", 
                              disable=not self.verbose,
                              file=sys.stdout)
            
            for chrom in chr_progress:
                chr_progress.set_description(f"Coverage analysis - {chrom}")
                
                # Multi-scale analysis
                for window_size in self.window_sizes:
                    candidates = self._analyze_chromosome_coverage(
                        bam, chrom, window_size
                    )
                    all_candidates.extend(candidates)
        
        return self._consolidate_candidates(all_candidates)
    
    def _analyze_chromosome_coverage(self, bam, chromosome, window_size):
        """Analyze coverage patterns for single chromosome at given window size"""
        chrom_length = bam.get_reference_length(chromosome)
        if chrom_length < window_size * 10:
            return []
        
        if self.verbose:
            print(f"  → Analyzing {chromosome} with {window_size}bp windows (length: {chrom_length:,}bp)")
        
        # Calculate coverage in windows
        coverage_data = self._calculate_windowed_coverage(
            bam, chromosome, window_size, chrom_length
        )
        
        if len(coverage_data) < 10:
            return []
        
        # Robust statistics using median
        coverage_values = np.array([w['coverage'] for w in coverage_data])
        positions = np.array([w['start'] for w in coverage_data])
        
        # Use median and MAD for robust threshold
        median_cov = np.median(coverage_values)
        mad_cov = median_abs_deviation(coverage_values)
        
        if median_cov < 1 or mad_cov == 0:
            return []
        
        # Dynamic threshold based on data distribution
        robust_threshold = median_cov + 3 * mad_cov
        fold_threshold = median_cov * self.min_fold_enrichment
        threshold = max(robust_threshold, fold_threshold, self.min_coverage)
        
        # Find peaks with local background subtraction
        candidates = self._find_coverage_peaks(
            coverage_values, positions, threshold, median_cov, window_size, chromosome
        )
        
        return candidates
    
    def _calculate_windowed_coverage(self, bam, chromosome, window_size, chrom_length):
        """Calculate coverage in sliding windows with local normalization"""
        coverage_data = []
        
        # Calculate number of windows for progress bar
        num_windows = (chrom_length + window_size - 1) // window_size
        
        # Create progress bar for window analysis
        window_progress = tqdm(
            range(0, chrom_length, window_size),
            desc=f"  → Processing {window_size}bp windows",
            total=num_windows,
            disable=not self.verbose,
            file=sys.stdout,
            leave=False
        )
        
        for start in window_progress:
            end = min(start + window_size, chrom_length)
            
            # Count reads with mapping quality filter
            count = 0
            for read in bam.fetch(chromosome, start, end):
                if read.mapping_quality >= 20 and not read.is_secondary:
                    count += 1
            
            # Normalize by window size
            normalized_coverage = count / (end - start) * window_size
            
            coverage_data.append({
                'start': start,
                'end': end,
                'coverage': normalized_coverage,
                'raw_count': count
            })
            
            # Update progress description periodically
            if len(coverage_data) % 1000 == 0:
                window_progress.set_postfix({
                    'windows': len(coverage_data),
                    'avg_cov': f"{np.mean([d['coverage'] for d in coverage_data[-100:]]):.1f}"
                })
        
        window_progress.close()
        return coverage_data
    
    def _find_coverage_peaks(self, coverage_values, positions, threshold, baseline, window_size, chromosome):
        """Find coverage peaks using local background subtraction"""
        candidates = []
        
        # Apply local background subtraction
        local_bg = self._calculate_local_background(coverage_values, window_size=10)
        adjusted_coverage = coverage_values - local_bg
        
        # Find peaks
        peak_indices, _ = find_peaks(
            adjusted_coverage,
            height=threshold - baseline,
            distance=max(5, 500 // window_size),
            width=max(2, 200 // window_size)
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
            
            # Filter by length
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
            
            # Calculate metrics
            fold_enrichment = region_median / max(flanking_bg, 1)
            uniformity = 1 - (region_mad / max(region_median, 1))
            
            # Create candidate WITHOUT confidence score (will be calculated by main detector)
            candidate = CircularCandidate(
                chromosome=chromosome,
                start=genomic_start,
                end=genomic_end,
                length=region_length,
                mean_coverage=region_median,
                fold_enrichment=fold_enrichment,
                coverage_uniformity=uniformity,
                confidence_score=0.0,  # Placeholder - will be set by ConfidenceScorer
                detection_method='coverage'
            )
            
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
                candidate.start <= current.end + 1000):  # Allow 1kb gap
                
                # Merge candidates - keep higher fold enrichment
                if candidate.fold_enrichment > current.fold_enrichment:
                    current = candidate
            else:
                consolidated.append(current)
                current = candidate
        
        consolidated.append(current)
        return consolidated