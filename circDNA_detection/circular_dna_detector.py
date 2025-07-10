#!/usr/bin/env python3
"""
ONT-Optimized Circular DNA Detection via Multi-Modal Analysis
Combines coverage patterns, junction detection, and split-read analysis
"""

import pysam
import numpy as np
import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Tuple
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation
import warnings
warnings.filterwarnings('ignore')

from .coverage_analyzer import CoverageAnalyzer
from .junction_detector import JunctionDetector
from .split_read_analyzer import SplitReadAnalyzer
from .utils import CircularCandidate, merge_overlapping_candidates, calculate_gc_content

class CircularDNADetector:
    def __init__(self, 
                 window_sizes=[100, 500, 1000],
                 min_fold_enrichment=1.5,
                 min_coverage=5,
                 min_length=200,
                 max_length=100000,
                 uniformity_threshold=0.4):
        """ONT-optimized circular DNA detector"""
        self.window_sizes = window_sizes
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.uniformity_threshold = uniformity_threshold
        
        # Initialize sub-analyzers
        self.coverage_analyzer = CoverageAnalyzer(
            window_sizes=window_sizes,
            min_fold_enrichment=min_fold_enrichment,
            min_coverage=min_coverage,
            uniformity_threshold=uniformity_threshold
        )
        self.junction_detector = JunctionDetector(min_support=3)
        self.split_read_analyzer = SplitReadAnalyzer(min_split_length=50)
    
    def detect_circular_dna(self, bam_file, reference_fasta, chromosome=None):
        """Multi-modal circular DNA detection pipeline"""
        print("Starting ONT-optimized circular DNA detection...")
        
        # 1. Coverage-based detection
        print("Phase 1: Coverage pattern analysis...")
        coverage_candidates = self.coverage_analyzer.detect_coverage_patterns(
            bam_file, chromosome
        )
        
        # 2. Junction detection
        print("Phase 2: Junction detection...")
        junction_candidates = self.junction_detector.detect_junctions(
            bam_file, chromosome
        )
        
        # 3. Split-read analysis
        print("Phase 3: Split-read analysis...")
        split_candidates = self.split_read_analyzer.analyze_split_reads(
            bam_file, chromosome
        )
        
        # 4. Candidate integration and scoring
        print("Phase 4: Candidate integration...")
        all_candidates = coverage_candidates + junction_candidates + split_candidates
        
        if not all_candidates:
            return []
        
        # Merge overlapping candidates
        merged_candidates = merge_overlapping_candidates(all_candidates)
        
        # Final scoring with multi-modal evidence
        final_candidates = self._score_integrated_candidates(
            merged_candidates, bam_file, reference_fasta
        )
        
        # Sort by confidence and apply final filters
        final_candidates = [c for c in final_candidates if c.confidence_score > 0.3]
        final_candidates.sort(key=lambda x: x.confidence_score, reverse=True)
        
        print(f"Detection complete: {len(final_candidates)} high-confidence candidates")
        return final_candidates
    
    def _score_integrated_candidates(self, candidates, bam_file, reference_fasta):
        """Score candidates using multi-modal evidence"""
        scored_candidates = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for candidate in candidates:
                # Multi-modal scoring
                scores = []
                
                # Coverage score (30%)
                if hasattr(candidate, 'fold_enrichment'):
                    coverage_score = min(1.0, candidate.fold_enrichment / 4.0)
                    scores.append(coverage_score * 0.3)
                
                # Junction score (25%)
                if hasattr(candidate, 'junction_support'):
                    junction_score = min(1.0, candidate.junction_support / 10.0)
                    scores.append(junction_score * 0.25)
                
                # Split-read score (20%)
                if hasattr(candidate, 'split_support'):
                    split_score = min(1.0, candidate.split_support / 5.0)
                    scores.append(split_score * 0.2)
                
                # Read orientation score (15%)
                orientation_score = self._calculate_orientation_score(
                    bam, candidate.chromosome, candidate.start, candidate.end
                )
                scores.append(orientation_score * 0.15)
                
                # Sequence complexity score (10%)
                if reference_fasta:
                    complexity_score = self._calculate_complexity_score(
                        reference_fasta, candidate.chromosome, candidate.start, candidate.end
                    )
                    scores.append(complexity_score * 0.1)
                
                candidate.confidence_score = sum(scores) if scores else 0.0
                
                if candidate.confidence_score > 0:
                    scored_candidates.append(candidate)
        
        return scored_candidates
    
    def _calculate_orientation_score(self, bam, chromosome, start, end):
        """Calculate read orientation bias score"""
        forward_reads = bam.count(chromosome, start, end, read_callback=lambda r: not r.is_reverse)
        reverse_reads = bam.count(chromosome, start, end, read_callback=lambda r: r.is_reverse)
        total_reads = forward_reads + reverse_reads
        
        if total_reads < 10:
            return 0.0
        
        # Circular DNA shows orientation bias
        ratio = min(forward_reads, reverse_reads) / total_reads
        return 1.0 - (2.0 * ratio)  # Higher score for more bias
    
    def _calculate_complexity_score(self, reference_fasta, chromosome, start, end):
        """Calculate sequence complexity score"""
        try:
            with pysam.FastaFile(reference_fasta) as fasta:
                sequence = fasta.fetch(chromosome, start, end)
                gc_content = calculate_gc_content(sequence)
                
                # Penalize extreme GC content
                gc_score = 1.0 - abs(gc_content - 0.5) * 2
                return max(0.0, gc_score)
        except:
            return 0.5
    
    def write_output(self, candidates, output_file):
        """Write results in enhanced BED format"""
        with open(output_file, 'w') as f:
            f.write("# ONT Circular DNA Detection Results\n")
            f.write("# chr\tstart\tend\tname\tscore\tstrand\tdetection_method\tconfidence\tdetails\n")
            
            for i, candidate in enumerate(candidates):
                name = f"circDNA_{i+1}"
                score = int(candidate.confidence_score * 1000)
                
                # Determine primary detection method
                methods = []
                if hasattr(candidate, 'fold_enrichment'):
                    methods.append(f"coverage({candidate.fold_enrichment:.1f})")
                if hasattr(candidate, 'junction_support'):
                    methods.append(f"junction({candidate.junction_support})")
                if hasattr(candidate, 'split_support'):
                    methods.append(f"split({candidate.split_support})")
                
                method_str = ";".join(methods)
                details = f"length={candidate.length};confidence={candidate.confidence_score:.3f}"
                
                f.write(f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
                       f"{name}\t{score}\t.\t{method_str}\t{candidate.confidence_score:.3f}\t{details}\n")

def main():
    parser = argparse.ArgumentParser(description="ONT-optimized circular DNA detection")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("reference_fasta", help="Reference FASTA file")
    parser.add_argument("-o", "--output", default="circular_dna_ont.bed", help="Output BED file")
    parser.add_argument("-c", "--chromosome", help="Analyze specific chromosome")
    parser.add_argument("--min-fold-enrichment", type=float, default=1.5)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--min-length", type=int, default=200)
    parser.add_argument("--max-length", type=int, default=100000)
    
    args = parser.parse_args()
    
    # Validate inputs
    try:
        with pysam.AlignmentFile(args.bam_file, "rb") as bam:
            if not bam.has_index():
                print("Warning: BAM file not indexed")
    except Exception as e:
        print(f"Error: Cannot read BAM file: {e}")
        sys.exit(1)
    
    # Initialize detector
    detector = CircularDNADetector(
        min_fold_enrichment=args.min_fold_enrichment,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        max_length=args.max_length
    )
    
    try:
        # Run detection
        candidates = detector.detect_circular_dna(
            args.bam_file, args.reference_fasta, args.chromosome
        )
        
        if candidates:
            detector.write_output(candidates, args.output)
            print(f"\nResults: {len(candidates)} candidates written to {args.output}")
            
            # Print top candidate
            top = candidates[0]
            print(f"Top candidate: {top.chromosome}:{top.start}-{top.end}")
            print(f"  Confidence: {top.confidence_score:.3f}")
            print(f"  Length: {top.length}bp")
        else:
            print("No circular DNA candidates detected")
            
    except Exception as e:
        print(f"Detection error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()