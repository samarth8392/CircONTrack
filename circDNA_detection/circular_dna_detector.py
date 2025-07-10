#!/usr/bin/env python3
"""
ONT-Optimized Circular DNA Detection via Multi-Modal Analysis
Combines coverage patterns, junction detection, and split-read analysis
WITH INTEGRATED CONFIDENCE SCORING SYSTEM
"""
import sys
import time
import logging
from tqdm import tqdm
import pysam
import numpy as np
from datetime import datetime
from typing import List, Dict, Tuple
from dataclasses import dataclass

# Import the confidence scoring components
from .confidence_scorer import ConfidenceScorer
from .coverage_analyzer import CoverageAnalyzer
from .utils import CircularCandidate

class CircularDNADetector:
    """Enhanced detector with comprehensive progress tracking, verbose output, and confidence scoring"""
    
    def __init__(self, min_fold_enrichment=1.5, min_coverage=5, min_length=200, 
                 max_length=100000, min_confidence=0.3, verbose=False, log_level='INFO'):
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.min_confidence = min_confidence
        self.verbose = verbose
        
        # Setup logging
        self.setup_logging(log_level)
        
        # Initialize confidence scorer and analyzers
        self.confidence_scorer = ConfidenceScorer()
        self.coverage_analyzer = CoverageAnalyzer(
            min_fold_enrichment=min_fold_enrichment,
            min_coverage=min_coverage
        )
        
        # Progress tracking
        self.total_steps = 0
        self.current_step = 0
        self.step_descriptions = []
        
    def setup_logging(self, log_level):
        """Setup logging configuration"""
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)
        
    def log_step(self, message, level='INFO'):
        """Log a step with progress information"""
        if self.verbose:
            progress = f"[{self.current_step}/{self.total_steps}]" if self.total_steps > 0 else ""
            getattr(self.logger, level.lower())(f"{progress} {message}")
            
    def detect_circular_dna(self, bam_file, reference_file, output_file=None, 
                          chromosome=None):
        """
        Main detection pipeline with comprehensive progress tracking and confidence scoring
        """
        start_time = time.time()
        self.log_step("Starting CircDNA detection pipeline with confidence scoring", "INFO")
        
        # Initialize progress tracking
        self._initialize_progress_tracking(bam_file, chromosome)
        
        try:
            # Phase 1: File validation and setup
            self._validate_inputs(bam_file, reference_file)
            
            # Phase 2: BAM file analysis
            bam_stats = self._analyze_bam_file(bam_file)
            
            # Phase 3: Run individual detection methods
            coverage_candidates = self._run_coverage_analysis(bam_file, chromosome)
            junction_candidates = self._run_junction_detection(bam_file, chromosome)
            split_candidates = self._run_split_read_analysis(bam_file, chromosome)
            
            # Phase 4: Multi-method integration with confidence scoring
            integrated_candidates = self._integrate_multi_method_results(
                coverage_candidates, junction_candidates, split_candidates
            )
            
            # Phase 5: Apply confidence filtering
            final_results = self._apply_confidence_filtering(integrated_candidates)
            
            # Phase 6: Final processing and output
            self._finalize_results(final_results, output_file)
            
            # Summary with confidence statistics
            self._print_confidence_summary(
                coverage_candidates, junction_candidates, split_candidates, 
                integrated_candidates, final_results, start_time
            )
            
            return final_results
            
        except Exception as e:
            self.logger.error(f"Error in detection pipeline: {str(e)}")
            raise
    
    def _initialize_progress_tracking(self, bam_file, chromosome):
        """Initialize progress tracking based on input"""
        self.log_step("Initializing progress tracking...")
        
        # Estimate total steps based on input
        if chromosome:
            self.total_steps = 12  # Single chromosome with confidence scoring
        else:
            # Estimate based on BAM file
            try:
                with pysam.AlignmentFile(bam_file, 'rb') as bam:
                    num_chromosomes = len(bam.references)
                    self.total_steps = 6 + (num_chromosomes * 2)  # Setup + per-chr steps
            except:
                self.total_steps = 60  # Default estimate
                
        self.current_step = 0
        self.log_step(f"Estimated {self.total_steps} total steps")
    
    def _validate_inputs(self, bam_file, reference_file):
        """Validate input files with progress updates"""
        self.current_step += 1
        self.log_step(f"Validating input files...")
        
        # Check BAM file
        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                if not bam.has_index():
                    self.log_step("WARNING: BAM file is not indexed, this will be slow", "WARNING")
                self.log_step(f"✓ BAM file validated: {bam_file}")
        except Exception as e:
            raise ValueError(f"Invalid BAM file: {str(e)}")
        
        # Check reference file
        try:
            with pysam.FastaFile(reference_file) as ref:
                num_contigs = len(ref.references)
                self.log_step(f"✓ Reference file validated: {reference_file} ({num_contigs} contigs)")
        except Exception as e:
            raise ValueError(f"Invalid reference file: {str(e)}")
    
    def _analyze_bam_file(self, bam_file):
        """Analyze BAM file statistics with progress"""
        self.current_step += 1
        self.log_step("Analyzing BAM file statistics...")
        
        stats = {
            'total_reads': 0,
            'mapped_reads': 0,
            'chromosomes': [],
            'avg_coverage': 0
        }
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            # Get chromosome information
            stats['chromosomes'] = list(bam.references)
            self.log_step(f"Found {len(stats['chromosomes'])} chromosomes/contigs")
            
            if self.verbose:
                # Quick read count (sample-based for speed)
                self.log_step("Sampling reads for statistics...")
                sample_size = 10000
                read_count = 0
                mapped_count = 0
                
                for i, read in enumerate(bam.fetch()):
                    if i >= sample_size:
                        break
                    read_count += 1
                    if not read.is_unmapped:
                        mapped_count += 1
                
                if read_count > 0:
                    stats['mapped_percentage'] = (mapped_count / read_count) * 100
                    self.log_step(f"Sample stats: {mapped_count}/{read_count} reads mapped "
                                f"({stats['mapped_percentage']:.1f}%)")
        
        return stats
    
    def _run_coverage_analysis(self, bam_file, chromosome):
        """Run coverage analysis with progress tracking"""
        self.current_step += 1
        self.log_step("Running coverage pattern analysis...")
        
        # Use the enhanced coverage analyzer
        candidates = self.coverage_analyzer.detect_coverage_patterns(bam_file, chromosome)
        
        # Calculate confidence scores for each candidate
        for candidate in candidates:
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        self.log_step(f"Coverage analysis complete: {len(candidates)} candidates found")
        return candidates
    
    def _run_junction_detection(self, bam_file, chromosome):
        """Run junction detection with progress tracking"""
        self.current_step += 1
        self.log_step("Running junction detection analysis...")
        
        candidates = []
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chr_name in chromosomes:
                self.log_step(f"  → Detecting junctions in {chr_name}")
                
                # Junction detection logic with confidence scoring
                junction_reads = []
                
                for read in bam.fetch(chr_name):
                    if self._is_junction_read(read):
                        junction_reads.append(read)
                
                # Group nearby junction reads
                junction_groups = self._group_junction_reads(junction_reads)
                
                # Create candidates from junction groups
                for group in junction_groups:
                    if len(group) >= 3:  # Minimum support
                        # Calculate junction boundaries
                        starts = [read.reference_start for read in group]
                        ends = [read.reference_end for read in group]
                        
                        candidate = CircularCandidate(
                            chromosome=chr_name,
                            start=min(starts),
                            end=max(ends),
                            length=max(ends) - min(starts),
                            junction_support=len(group),
                            confidence_score=0.0,
                            detection_method='junction'
                        )
                        
                        # Calculate confidence score
                        candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
                        
                        # Apply length filters
                        if self.min_length <= candidate.length <= self.max_length:
                            candidates.append(candidate)
        
        self.log_step(f"Junction detection complete: {len(candidates)} candidates found")
        return candidates
    
    def _run_split_read_analysis(self, bam_file, chromosome):
        """Run split-read analysis with progress tracking"""
        self.current_step += 1
        self.log_step("Running split-read analysis...")
        
        candidates = []
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chr_name in chromosomes:
                self.log_step(f"  → Analyzing split reads in {chr_name}")
                
                # Split-read detection logic
                split_reads = []
                
                for read in bam.fetch(chr_name):
                    if self._is_split_read(read):
                        split_reads.append(read)
                
                # Group nearby split reads
                split_groups = self._group_split_reads(split_reads)
                
                # Create candidates from split groups
                for group in split_groups:
                    if len(group) >= 3:  # Minimum support
                        # Calculate split boundaries
                        starts = [read.reference_start for read in group]
                        ends = [read.reference_end for read in group]
                        
                        candidate = CircularCandidate(
                            chromosome=chr_name,
                            start=min(starts),
                            end=max(ends),
                            length=max(ends) - min(starts),
                            split_support=len(group),
                            confidence_score=0.0,
                            detection_method='split_read'
                        )
                        
                        # Calculate confidence score
                        candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
                        
                        # Apply length filters
                        if self.min_length <= candidate.length <= self.max_length:
                            candidates.append(candidate)
        
        self.log_step(f"Split-read analysis complete: {len(candidates)} candidates found")
        return candidates
    
    def _integrate_multi_method_results(self, coverage_candidates, junction_candidates, split_candidates):
        """Integrate results from multiple detection methods"""
        self.current_step += 1
        self.log_step("Integrating multi-method results...")
        
        all_candidates = coverage_candidates + junction_candidates + split_candidates
        
        if not all_candidates:
            return []
        
        # Sort candidates by chromosome and position
        all_candidates.sort(key=lambda x: (x.chromosome, x.start))
        
        # Merge overlapping candidates
        merged_candidates = []
        current = all_candidates[0]
        
        for candidate in all_candidates[1:]:
            if self._should_merge_candidates(current, candidate):
                current = self._merge_two_candidates(current, candidate)
            else:
                merged_candidates.append(current)
                current = candidate
        
        merged_candidates.append(current)
        
        # Recalculate confidence scores for merged candidates
        for candidate in merged_candidates:
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
        
        self.log_step(f"Multi-method integration complete: {len(merged_candidates)} integrated candidates")
        return merged_candidates
    
    def _apply_confidence_filtering(self, candidates):
        """Apply confidence score filtering"""
        self.current_step += 1
        self.log_step(f"Applying confidence filtering (threshold: {self.min_confidence:.2f})...")
        
        # Filter by confidence threshold
        high_confidence_candidates = [
            c for c in candidates if c.confidence_score >= self.min_confidence
        ]
        
        # Sort by confidence score (descending)
        high_confidence_candidates.sort(key=lambda x: x.confidence_score, reverse=True)
        
        self.log_step(f"Confidence filtering complete: {len(high_confidence_candidates)} high-confidence candidates")
        return high_confidence_candidates
    
    def _finalize_results(self, candidates, output_file):
        """Finalize results and write output"""
        self.current_step += 1
        self.log_step("Finalizing results...")
        
        if output_file:
            self._write_output_with_confidence(candidates, output_file)
        
        self.log_step(f"Results finalized: {len(candidates)} total candidates")
        return candidates
    
    def _write_output_with_confidence(self, candidates, output_file):
        """Write results to BED file with confidence scores"""
        self.log_step(f"Writing results to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("# CircDNA Detection Results with Confidence Scoring\n")
            f.write("# chr\tstart\tend\tname\tconfidence\tstrand\tmethod\tlength\tevidence\n")
            
            for i, candidate in enumerate(candidates):
                # Build evidence string
                evidence_parts = []
                if candidate.fold_enrichment:
                    evidence_parts.append(f"fold={candidate.fold_enrichment:.2f}")
                if candidate.mean_coverage:
                    evidence_parts.append(f"cov={candidate.mean_coverage:.1f}")
                if candidate.junction_support:
                    evidence_parts.append(f"junc={candidate.junction_support}")
                if candidate.split_support:
                    evidence_parts.append(f"split={candidate.split_support}")
                
                evidence_str = ",".join(evidence_parts) if evidence_parts else "."
                
                f.write(f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
                       f"circDNA_{i+1}\t{candidate.confidence_score:.3f}\t.\t"
                       f"{candidate.detection_method}\t{candidate.length}\t{evidence_str}\n")
    
    def _print_confidence_summary(self, coverage_cands, junction_cands, split_cands, 
                                integrated_cands, final_cands, start_time):
        """Print comprehensive confidence scoring summary"""
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*70)
        print("CIRCULAR DNA DETECTION SUMMARY WITH CONFIDENCE SCORING")
        print("="*70)
        
        print(f"Analysis completed in: {elapsed_time:.2f} seconds")
        print(f"Confidence threshold: {self.min_confidence:.2f}")
        
        print(f"\nRaw candidates found:")
        print(f"  Coverage-based: {len(coverage_cands)}")
        print(f"  Junction-based: {len(junction_cands)}")
        print(f"  Split-read based: {len(split_cands)}")
        print(f"  Total raw: {len(coverage_cands) + len(junction_cands) + len(split_cands)}")
        
        print(f"\nAfter multi-method integration: {len(integrated_cands)}")
        print(f"After confidence filtering: {len(final_cands)}")
        
        if final_cands:
            # Method breakdown
            method_counts = {}
            for result in final_cands:
                method = result.detection_method
                method_counts[method] = method_counts.get(method, 0) + 1
            
            print(f"\nHigh-confidence candidates by detection method:")
            for method, count in method_counts.items():
                print(f"  {method}: {count}")
            
            # Confidence distribution
            confidences = [c.confidence_score for c in final_cands]
            print(f"\nConfidence score distribution:")
            print(f"  Mean: {np.mean(confidences):.3f}")
            print(f"  Median: {np.median(confidences):.3f}")
            print(f"  Range: {min(confidences):.3f} - {max(confidences):.3f}")
            
            # Top candidates
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
                if candidate.mean_coverage:
                    evidence.append(f"coverage={candidate.mean_coverage:.1f}")
                if candidate.junction_support:
                    evidence.append(f"junctions={candidate.junction_support}")
                if candidate.split_support:
                    evidence.append(f"splits={candidate.split_support}")
                
                if evidence:
                    print(f"     Evidence: {', '.join(evidence)}")
                print()
        
        print("="*70)
    
    def _should_merge_candidates(self, c1, c2):
        """Check if two candidates should be merged"""
        if c1.chromosome != c2.chromosome:
            return False
        
        # Allow 1kb gap for merging
        overlap_distance = 1000
        overlap = (c2.start <= c1.end + overlap_distance and
                  c1.start <= c2.end + overlap_distance)
        
        return overlap
    
    def _merge_two_candidates(self, c1, c2):
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
    
    def _group_junction_reads(self, junction_reads):
        """Group nearby junction reads"""
        if not junction_reads:
            return []
        
        # Sort by position
        junction_reads.sort(key=lambda r: r.reference_start)
        
        groups = []
        current_group = [junction_reads[0]]
        
        for read in junction_reads[1:]:
            # Group reads within 1kb
            if read.reference_start - current_group[-1].reference_start <= 1000:
                current_group.append(read)
            else:
                groups.append(current_group)
                current_group = [read]
        
        groups.append(current_group)
        return groups
    
    def _group_split_reads(self, split_reads):
        """Group nearby split reads"""
        if not split_reads:
            return []
        
        # Sort by position
        split_reads.sort(key=lambda r: r.reference_start)
        
        groups = []
        current_group = [split_reads[0]]
        
        for read in split_reads[1:]:
            # Group reads within 1kb
            if read.reference_start - current_group[-1].reference_start <= 1000:
                current_group.append(read)
            else:
                groups.append(current_group)
                current_group = [read]
        
        groups.append(current_group)
        return groups
    
    def _is_junction_read(self, read):
        """Check if read shows junction signature"""
        # Enhanced junction detection logic
        if read.is_unmapped or read.mapping_quality < 20:
            return False
        
        # Check for soft clipping indicating junction
        if read.cigarstring and 'S' in read.cigarstring:
            # Parse CIGAR to check for significant soft clipping
            import re
            soft_clips = re.findall(r'(\d+)S', read.cigarstring)
            if soft_clips:
                max_clip = max(int(clip) for clip in soft_clips)
                return max_clip >= 20  # Minimum 20bp soft clip
        
        return False
    
    def _is_split_read(self, read):
        """Check if read is split alignment"""
        # Enhanced split-read detection
        if read.is_unmapped or read.mapping_quality < 20:
            return False
        
        # Check for supplementary alignment tag
        if read.has_tag('SA'):
            return True
        
        # Check for multiple large soft clips
        if read.cigarstring and read.cigarstring.count('S') >= 2:
            import re
            soft_clips = re.findall(r'(\d+)S', read.cigarstring)
            if len(soft_clips) >= 2:
                return any(int(clip) >= 20 for clip in soft_clips)
        
        return False


def main():
    """Enhanced main function with confidence scoring parameters"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='CircDNA Detection with Confidence Scoring System'
    )
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('reference_file', help='Reference FASTA file')
    parser.add_argument('-o', '--output', help='Output BED file', 
                       default='circular_dna_confident.bed')
    parser.add_argument('-c', '--chromosome', help='Analyze specific chromosome only')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Enable verbose output')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                       default='INFO', help='Set logging level')
    
    # Detection parameters
    parser.add_argument('--min-fold-enrichment', type=float, default=1.5,
                       help='Minimum fold enrichment for coverage peaks')
    parser.add_argument('--min-coverage', type=int, default=5,
                       help='Minimum coverage threshold')
    parser.add_argument('--min-length', type=int, default=200,
                       help='Minimum candidate length')
    parser.add_argument('--max-length', type=int, default=100000,
                       help='Maximum candidate length')
    
    # Confidence scoring parameters
    parser.add_argument('--min-confidence', type=float, default=0.3,
                       help='Minimum confidence score threshold')

    args = parser.parse_args()
    
    # Create detector with confidence scoring
    detector = CircularDNADetector(
        min_fold_enrichment=args.min_fold_enrichment,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        max_length=args.max_length,
        min_confidence=args.min_confidence,
        verbose=args.verbose,
        log_level=args.log_level
    )
    
    # Run detection
    results = detector.detect_circular_dna(
        args.bam_file,
        args.reference_file,
        args.output,
        args.chromosome
    )
    
    print(f"\nAnalysis complete! Results written to {args.output}")
    print(f"Found {len(results)} high-confidence circular DNA candidates")


if __name__ == "__main__":
    main()