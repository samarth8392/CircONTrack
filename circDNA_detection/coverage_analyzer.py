#!/usr/bin/env python3
"""
Optimized Coverage Pattern Analysis
Streamlined for integration with main CircDNA detector
"""

import pysam
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation
from typing import List, Dict, Tuple, Optional
import logging


class CoverageAnalyzer:
    """
    Optimized coverage pattern analyzer with improved memory efficiency
    and better integration with the main detection pipeline
    """
    
    def __init__(self, window_sizes=[500, 1000, 2000], 
                 min_fold_enrichment=1.5, min_coverage=5, 
                 uniformity_threshold=0.4, mapping_quality_threshold=20):
        """
        Initialize coverage analyzer with configurable parameters
        
        Args:
            window_sizes: List of window sizes for multi-scale analysis
            min_fold_enrichment: Minimum fold enrichment over background
            min_coverage: Minimum absolute coverage threshold
            uniformity_threshold: Minimum coverage uniformity score
            mapping_quality_threshold: Minimum mapping quality for reads
        """
        self.window_sizes = window_sizes
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.uniformity_threshold = uniformity_threshold
        self.mapping_quality_threshold = mapping_quality_threshold
        
        # Setup logging
        self.logger = logging.getLogger(__name__)
    
    def analyze_chromosome_coverage(self, bam: pysam.AlignmentFile, 
                                   chromosome: str, chr_length: int) -> List[Dict]:
        """
        Analyze coverage patterns for a single chromosome using multi-scale approach
        
        Args:
            bam: Opened BAM file handle
            chromosome: Chromosome name
            chr_length: Length of chromosome
            
        Returns:
            List of coverage-based candidate dictionaries
        """
        all_candidates = []
        
        # Skip very short chromosomes
        if chr_length < max(self.window_sizes) * 10:
            self.logger.warning(f"Chromosome {chromosome} too short for analysis")
            return []
        
        # Multi-scale analysis
        for window_size in self.window_sizes:
            try:
                candidates = self._analyze_single_scale(
                    bam, chromosome, chr_length, window_size
                )
                all_candidates.extend(candidates)
                
                self.logger.debug(f"Found {len(candidates)} candidates at {window_size}bp scale")
                
            except Exception as e:
                self.logger.error(f"Error analyzing {chromosome} at {window_size}bp: {str(e)}")
                continue
        
        # Consolidate overlapping candidates from different scales
        return self._consolidate_multi_scale_candidates(all_candidates)
    
    def _analyze_single_scale(self, bam: pysam.AlignmentFile, 
                             chromosome: str, chr_length: int, 
                             window_size: int) -> List[Dict]:
        """
        Analyze coverage at a single window scale
        
        Args:
            bam: BAM file handle
            chromosome: Chromosome name
            chr_length: Chromosome length
            window_size: Window size for analysis
            
        Returns:
            List of candidate dictionaries
        """
        # Calculate windowed coverage
        coverage_data = self._calculate_windowed_coverage(
            bam, chromosome, chr_length, window_size
        )
        
        if len(coverage_data) < 10:
            return []
        
        # Extract coverage values and positions
        coverage_values = np.array([w['coverage'] for w in coverage_data])
        positions = np.array([w['start'] for w in coverage_data])
        
        # Calculate robust statistics
        baseline_stats = self._calculate_baseline_statistics(coverage_values)
        
        if baseline_stats['median'] < 1:
            return []
        
        # Find enriched regions
        candidates = self._find_enriched_regions(
            coverage_values, positions, baseline_stats, 
            window_size, chromosome
        )
        
        return candidates
    
    def _calculate_windowed_coverage(self, bam: pysam.AlignmentFile, 
                                   chromosome: str, chr_length: int, 
                                   window_size: int) -> List[Dict]:
        """
        Calculate coverage in non-overlapping windows with quality filtering
        
        Args:
            bam: BAM file handle
            chromosome: Chromosome name
            chr_length: Chromosome length
            window_size: Size of analysis windows
            
        Returns:
            List of window coverage data
        """
        coverage_data = []
        
        for start in range(0, chr_length, window_size):
            end = min(start + window_size, chr_length)
            window_length = end - start
            
            # Count high-quality reads in window
            read_count = 0
            total_bases = 0
            
            try:
                for read in bam.fetch(chromosome, start, end):
                    if self._is_valid_read(read):
                        # Calculate overlap with window
                        overlap_start = max(read.reference_start, start)
                        overlap_end = min(read.reference_end, end)
                        overlap_length = max(0, overlap_end - overlap_start)
                        
                        read_count += 1
                        total_bases += overlap_length
                
                # Calculate normalized coverage
                coverage = total_bases / window_length if window_length > 0 else 0
                
                coverage_data.append({
                    'start': start,
                    'end': end,
                    'length': window_length,
                    'coverage': coverage,
                    'read_count': read_count
                })
                
            except Exception as e:
                self.logger.warning(f"Error processing window {start}-{end}: {str(e)}")
                continue
        
        return coverage_data
    
    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """
        Check if read passes quality filters
        
        Args:
            read: BAM read alignment
            
        Returns:
            True if read passes filters
        """
        return (not read.is_unmapped and
                not read.is_secondary and
                not read.is_supplementary and
                read.mapping_quality >= self.mapping_quality_threshold)
    
    def _calculate_baseline_statistics(self, coverage_values: np.ndarray) -> Dict:
        """
        Calculate robust baseline statistics using median and MAD
        
        Args:
            coverage_values: Array of coverage values
            
        Returns:
            Dictionary with baseline statistics
        """
        median_cov = np.median(coverage_values)
        mad_cov = median_abs_deviation(coverage_values)
        
        # Calculate various thresholds
        robust_threshold = median_cov + 3 * mad_cov
        fold_threshold = median_cov * self.min_fold_enrichment
        final_threshold = max(robust_threshold, fold_threshold, self.min_coverage)
        
        return {
            'median': median_cov,
            'mad': mad_cov,
            'robust_threshold': robust_threshold,
            'fold_threshold': fold_threshold,
            'final_threshold': final_threshold
        }
    
    def _find_enriched_regions(self, coverage_values: np.ndarray, 
                              positions: np.ndarray, baseline_stats: Dict,
                              window_size: int, chromosome: str) -> List[Dict]:
        """
        Find enriched regions using peak detection with local background correction
        
        Args:
            coverage_values: Array of coverage values
            positions: Array of genomic positions
            baseline_stats: Baseline statistics dictionary
            window_size: Window size used
            chromosome: Chromosome name
            
        Returns:
            List of candidate dictionaries
        """
        candidates = []
        
        # Apply local background subtraction
        local_background = self._calculate_local_background(coverage_values)
        adjusted_coverage = coverage_values - local_background
        
        # Find peaks using scipy
        min_peak_distance = max(3, 1000 // window_size)  # Minimum 1kb separation
        min_peak_width = max(1, 200 // window_size)      # Minimum 200bp width
        
        peak_indices, peak_properties = find_peaks(
            adjusted_coverage,
            height=baseline_stats['final_threshold'] - baseline_stats['median'],
            distance=min_peak_distance,
            width=min_peak_width,
            prominence=baseline_stats['mad']
        )
        
        # Process each peak
        for i, peak_idx in enumerate(peak_indices):
            try:
                candidate = self._create_peak_candidate(
                    peak_idx, coverage_values, positions, 
                    baseline_stats, window_size, chromosome,
                    peak_properties, i
                )
                
                if candidate:
                    candidates.append(candidate)
                    
            except Exception as e:
                self.logger.warning(f"Error processing peak {i}: {str(e)}")
                continue
        
        return candidates
    
    def _calculate_local_background(self, coverage_values: np.ndarray, 
                                  background_window: int = 20) -> np.ndarray:
        """
        Calculate local background using rolling median
        
        Args:
            coverage_values: Array of coverage values
            background_window: Size of background calculation window
            
        Returns:
            Array of local background values
        """
        local_bg = np.zeros_like(coverage_values)
        half_window = background_window // 2
        
        for i in range(len(coverage_values)):
            start_idx = max(0, i - half_window)
            end_idx = min(len(coverage_values), i + half_window + 1)
            local_bg[i] = np.median(coverage_values[start_idx:end_idx])
        
        return local_bg
    
    def _create_peak_candidate(self, peak_idx: int, coverage_values: np.ndarray,
                              positions: np.ndarray, baseline_stats: Dict,
                              window_size: int, chromosome: str,
                              peak_properties: Dict, peak_id: int) -> Optional[Dict]:
        """
        Create a candidate dictionary from peak information
        
        Args:
            peak_idx: Index of peak in coverage array
            coverage_values: Array of coverage values
            positions: Array of genomic positions
            baseline_stats: Baseline statistics
            window_size: Window size used
            chromosome: Chromosome name
            peak_properties: Properties from peak detection
            peak_id: Unique peak identifier
            
        Returns:
            Candidate dictionary or None if invalid
        """
        # Find peak boundaries
        start_idx, end_idx = self._find_peak_boundaries(
            coverage_values, peak_idx, baseline_stats['median']
        )
        
        # Convert to genomic coordinates
        genomic_start = positions[start_idx]
        genomic_end = positions[min(end_idx, len(positions) - 1)] + window_size
        region_length = genomic_end - genomic_start
        
        # Filter by length constraints
        if region_length < 200 or region_length > 100000:
            return None
        
        # Calculate region statistics
        region_coverage = coverage_values[start_idx:end_idx + 1]
        region_stats = self._calculate_region_statistics(
            region_coverage, coverage_values, start_idx, end_idx, baseline_stats
        )
        
        # Create candidate dictionary
        candidate = {
            'chromosome': chromosome,
            'start': genomic_start,
            'end': genomic_end,
            'length': region_length,
            'mean_coverage': region_stats['mean_coverage'],
            'fold_enrichment': region_stats['fold_enrichment'],
            'coverage_uniformity': region_stats['uniformity'],
            'detection_method': 'coverage',
            'window_size': window_size,
            'peak_height': coverage_values[peak_idx],
            'peak_id': f"{chromosome}_{peak_id}_{window_size}"
        }
        
        return candidate
    
    def _find_peak_boundaries(self, coverage_values: np.ndarray, 
                             peak_idx: int, threshold: float) -> Tuple[int, int]:
        """
        Find peak boundaries using threshold crossing
        
        Args:
            coverage_values: Array of coverage values
            peak_idx: Index of peak center
            threshold: Threshold for boundary detection
            
        Returns:
            Tuple of (start_index, end_index)
        """
        # Find left boundary
        left_boundary = peak_idx
        for i in range(peak_idx - 1, -1, -1):
            if coverage_values[i] <= threshold:
                left_boundary = i + 1
                break
            left_boundary = i
        
        # Find right boundary
        right_boundary = peak_idx
        for i in range(peak_idx + 1, len(coverage_values)):
            if coverage_values[i] <= threshold:
                right_boundary = i - 1
                break
            right_boundary = i
        
        return left_boundary, right_boundary
    
    def _calculate_region_statistics(self, region_coverage: np.ndarray,
                                   all_coverage: np.ndarray,
                                   start_idx: int, end_idx: int,
                                   baseline_stats: Dict) -> Dict:
        """
        Calculate comprehensive statistics for a candidate region
        
        Args:
            region_coverage: Coverage values in the region
            all_coverage: All coverage values
            start_idx: Start index of region
            end_idx: End index of region
            baseline_stats: Baseline statistics
            
        Returns:
            Dictionary with region statistics
        """
        # Basic statistics
        mean_coverage = np.mean(region_coverage)
        median_coverage = np.median(region_coverage)
        std_coverage = np.std(region_coverage)
        
        # Calculate flanking background
        flanking_bg = self._calculate_flanking_background(
            all_coverage, start_idx, end_idx
        )
        
        # Calculate fold enrichment
        fold_enrichment = mean_coverage / max(flanking_bg, baseline_stats['median'], 1)
        
        # Calculate uniformity (1 - coefficient of variation)
        uniformity = 1 - (std_coverage / max(mean_coverage, 1)) if mean_coverage > 0 else 0
        uniformity = max(0, uniformity)  # Ensure non-negative
        
        return {
            'mean_coverage': mean_coverage,
            'median_coverage': median_coverage,
            'std_coverage': std_coverage,
            'fold_enrichment': fold_enrichment,
            'uniformity': uniformity,
            'flanking_background': flanking_bg
        }
    
    def _calculate_flanking_background(self, coverage_values: np.ndarray,
                                     start_idx: int, end_idx: int) -> float:
        """
        Calculate background coverage from flanking regions
        
        Args:
            coverage_values: All coverage values
            start_idx: Start index of region
            end_idx: End index of region
            
        Returns:
            Flanking background coverage
        """
        region_size = end_idx - start_idx
        flank_size = min(region_size, 30)  # Adaptive flank size
        
        # Left flank
        left_start = max(0, start_idx - flank_size)
        left_flank = coverage_values[left_start:start_idx]
        
        # Right flank
        right_end = min(len(coverage_values), end_idx + flank_size)
        right_flank = coverage_values[end_idx:right_end]
        
        # Combine and calculate median
        flanking_values = np.concatenate([left_flank, right_flank])
        
        return np.median(flanking_values) if len(flanking_values) > 0 else 0
    
    def _consolidate_multi_scale_candidates(self, candidates: List[Dict]) -> List[Dict]:
        """
        Consolidate candidates from multiple scales, keeping the best representatives
        
        Args:
            candidates: List of candidate dictionaries from all scales
            
        Returns:
            List of consolidated candidates
        """
        if not candidates:
            return []
        
        # Sort by position for efficient processing
        candidates.sort(key=lambda x: (x['chromosome'], x['start']))
        
        consolidated = []
        current_group = [candidates[0]]
        
        for candidate in candidates[1:]:
            # Check if candidate overlaps with current group
            if self._overlaps_with_group(candidate, current_group):
                current_group.append(candidate)
            else:
                # Process current group and start new one
                best_candidate = self._select_best_candidate(current_group)
                consolidated.append(best_candidate)
                current_group = [candidate]
        
        # Process final group
        if current_group:
            best_candidate = self._select_best_candidate(current_group)
            consolidated.append(best_candidate)
        
        return consolidated
    
    def _overlaps_with_group(self, candidate: Dict, group: List[Dict]) -> bool:
        """
        Check if candidate overlaps with any candidate in the group
        
        Args:
            candidate: Candidate to check
            group: List of candidates in current group
            
        Returns:
            True if candidate overlaps with group
        """
        for existing in group:
            if (candidate['chromosome'] == existing['chromosome'] and
                not (candidate['end'] < existing['start'] or 
                     candidate['start'] > existing['end'])):
                return True
        return False
    
    def _select_best_candidate(self, group: List[Dict]) -> Dict:
        """
        Select the best candidate from overlapping group
        
        Args:
            group: List of overlapping candidates
            
        Returns:
            Best candidate from the group
        """
        if len(group) == 1:
            return group[0]
        
        # Score candidates based on multiple criteria
        best_candidate = None
        best_score = -1
        
        for candidate in group:
            # Multi-criteria scoring
            fold_score = min(candidate['fold_enrichment'] / 5.0, 1.0)
            uniformity_score = candidate['coverage_uniformity']
            coverage_score = min(candidate['mean_coverage'] / 20.0, 1.0)
            
            # Combined score with weights
            combined_score = (fold_score * 0.4 + 
                            uniformity_score * 0.3 + 
                            coverage_score * 0.3)
            
            if combined_score > best_score:
                best_score = combined_score
                best_candidate = candidate
        
        return best_candidate if best_candidate else group[0]