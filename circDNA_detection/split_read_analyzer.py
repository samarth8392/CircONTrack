#!/usr/bin/env python3
"""
ONT-Optimized Split Read Analysis for Circular DNA Detection
Handles noise and characteristics specific to long-read sequencing
"""

import pysam
import numpy as np
from collections import defaultdict
import re
from .utils import CircularCandidate
from .confidence_scorer import ConfidenceScorer

class SplitReadAnalyzer:
    def __init__(self, min_split_length=100, min_support=3, max_distance=2000,
                 min_alignment_quality=30, min_split_mapq=20,
                 max_nm_rate=0.15, min_overlap_ratio=0.8, verbose=False):
        """
        Initialize split read analyzer for ONT data
        
        Parameters:
        -----------
        min_split_length : int
            Minimum length for split alignment (increased for ONT)
        min_support : int
            Minimum number of supporting reads
        max_distance : int
            Maximum distance between splits to consider circular
        min_alignment_quality : int
            Minimum alignment score for split alignments
        min_split_mapq : int
            Minimum mapping quality for split alignments
        max_nm_rate : float
            Maximum edit distance rate (NM/length)
        min_overlap_ratio : float
            Minimum ratio of combined alignment length to read length
        verbose : bool
            Enable verbose output
        """
        self.min_split_length = min_split_length
        self.min_support = min_support
        self.max_distance = max_distance
        self.min_alignment_quality = min_alignment_quality
        self.min_split_mapq = min_split_mapq
        self.max_nm_rate = max_nm_rate
        self.min_overlap_ratio = min_overlap_ratio
        self.verbose = verbose
        self.confidence_scorer = ConfidenceScorer()
    
    def analyze_split_reads(self, bam_file, chromosome=None):
        """Analyze split reads to detect circular DNA junctions"""
        print("  Analyzing split reads...")
        
        split_candidates = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chrom in chromosomes:
                # Collect high-quality split alignments
                splits = self._collect_split_alignments(bam, chrom)
                
                # Filter noise and validate patterns
                validated_splits = self._validate_split_patterns(splits, chrom)
                
                # Create candidates from validated splits
                candidates = self._create_split_candidates(validated_splits, chrom)
                split_candidates.extend(candidates)
        
        return split_candidates
    
    def _collect_split_alignments(self, bam, chromosome):
        """Collect split read alignments with strict quality filtering"""
        split_reads = defaultdict(list)
        
        for read in bam.fetch(chromosome):
            if read.is_unmapped or read.is_secondary:
                continue
            
            # Primary alignment quality check
            if read.mapping_quality < self.min_split_mapq:
                continue
            
            # Check for supplementary alignments
            if not read.has_tag('SA'):
                continue
            
            # Get read length for validation
            read_length = read.infer_read_length()
            if not read_length:
                continue
            
            # Parse SA tag
            sa_tag = read.get_tag('SA')
            supplementary_alns = self._parse_sa_tag(sa_tag)
            
            # Create primary alignment info
            primary_aln = {
                'chromosome': chromosome,
                'start': read.reference_start,
                'end': read.reference_end,
                'strand': '-' if read.is_reverse else '+',
                'mapq': read.mapping_quality,
                'length': read.reference_end - read.reference_start,
                'nm': read.get_tag('NM') if read.has_tag('NM') else 0,
                'cigar': read.cigarstring
            }
            
            # Validate and filter supplementary alignments
            valid_supplements = []
            for supp in supplementary_alns:
                if self._is_valid_supplementary(supp, primary_aln, read_length):
                    valid_supplements.append(supp)
            
            if valid_supplements:
                split_reads[read.query_name] = {
                    'primary': primary_aln,
                    'supplementary': valid_supplements,
                    'read_length': read_length,
                    'read_name': read.query_name
                }
        
        return split_reads
    
    def _parse_sa_tag(self, sa_tag):
        """Parse SA tag with error handling"""
        supplements = []
        
        for alignment in sa_tag.strip().split(';'):
            if not alignment:
                continue
            
            parts = alignment.split(',')
            if len(parts) >= 6:
                try:
                    chromosome = parts[0]
                    pos = int(parts[1]) - 1  # Convert to 0-based
                    strand = parts[2]
                    cigar = parts[3]
                    mapq = int(parts[4])
                    nm = int(parts[5]) if parts[5].isdigit() else 0
                    
                    # Calculate alignment length and stats
                    aln_stats = self._parse_cigar_stats(cigar)
                    
                    if aln_stats['ref_length'] >= self.min_split_length:
                        supplements.append({
                            'chromosome': chromosome,
                            'start': pos,
                            'end': pos + aln_stats['ref_length'],
                            'strand': strand,
                            'cigar': cigar,
                            'mapq': mapq,
                            'nm': nm,
                            'length': aln_stats['ref_length'],
                            'aligned_length': aln_stats['aligned_length'],
                            'clipped': aln_stats['clipped']
                        })
                except (ValueError, IndexError):
                    continue
        
        return supplements
    
    def _parse_cigar_stats(self, cigar_string):
        """Parse CIGAR string to get alignment statistics"""
        ref_length = 0
        aligned_length = 0
        clipped = 0
        
        matches = re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)
        
        for count, op in matches:
            count = int(count)
            if op in 'M=X':  # Alignment match/mismatch
                ref_length += count
                aligned_length += count
            elif op in 'DN':  # Deletion/Reference skip
                ref_length += count
            elif op == 'I':  # Insertion
                aligned_length += count
            elif op in 'SH':  # Soft/Hard clip
                clipped += count
        
        return {
            'ref_length': ref_length,
            'aligned_length': aligned_length,
            'clipped': clipped
        }
    
    def _is_valid_supplementary(self, supp, primary, read_length):
        """Validate supplementary alignment quality"""
        # Basic filters
        if supp['chromosome'] != primary['chromosome']:
            return False
        
        if supp['mapq'] < self.min_split_mapq:
            return False
        
        if supp['length'] < self.min_split_length:
            return False
        
        # Edit distance rate check
        if supp['nm'] > 0 and supp['aligned_length'] > 0:
            nm_rate = supp['nm'] / supp['aligned_length']
            if nm_rate > self.max_nm_rate:
                return False
        
        # Check for sufficient coverage of read
        total_aligned = primary['length'] + supp['length']
        coverage_ratio = total_aligned / read_length
        
        if coverage_ratio < self.min_overlap_ratio:
            return False
        
        # Check for overlap between primary and supplementary
        # (overlaps often indicate noise rather than true splits)
        if (primary['start'] < supp['end'] and supp['start'] < primary['end']):
            overlap = min(primary['end'], supp['end']) - max(primary['start'], supp['start'])
            if overlap > 50:  # Significant overlap
                return False
        
        return True
    
    def _validate_split_patterns(self, split_reads, chromosome):
        """Validate split patterns for circular DNA signatures"""
        validated_splits = []
        
        for read_name, split_info in split_reads.items():
            primary = split_info['primary']
            
            for supp in split_info['supplementary']:
                # Check for circular DNA patterns
                pattern = self._check_circular_pattern(primary, supp)
                
                if pattern:
                    validated_splits.append({
                        'read_name': read_name,
                        'pattern': pattern,
                        'primary': primary,
                        'supplementary': supp,
                        'quality_score': self._calculate_split_quality(primary, supp)
                    })
        
        return validated_splits
    
    def _check_circular_pattern(self, primary, supp):
        """Check if split alignment suggests circular DNA"""
        # Pattern 1: Back-to-back alignment (suggesting circular junction)
        # Primary end is close to supplementary start
        if abs(primary['end'] - supp['start']) <= self.max_distance:
            if primary['strand'] == supp['strand']:
                return {
                    'type': 'back_to_back',
                    'junction_pos': (primary['end'] + supp['start']) // 2,
                    'distance': abs(primary['end'] - supp['start']),
                    'orientation': 'consistent'
                }
        
        # Pattern 2: Primary start close to supplementary end
        if abs(primary['start'] - supp['end']) <= self.max_distance:
            if primary['strand'] == supp['strand']:
                return {
                    'type': 'back_to_back',
                    'junction_pos': (primary['start'] + supp['end']) // 2,
                    'distance': abs(primary['start'] - supp['end']),
                    'orientation': 'consistent'
                }
        
        # Pattern 3: Inverted repeats (common in circular DNA)
        if primary['strand'] != supp['strand']:
            # Check if ends are close (inverted duplication)
            if abs(primary['end'] - supp['end']) <= self.max_distance:
                return {
                    'type': 'inverted',
                    'junction_pos': (primary['end'] + supp['end']) // 2,
                    'distance': abs(primary['end'] - supp['end']),
                    'orientation': 'inverted'
                }
            
            if abs(primary['start'] - supp['start']) <= self.max_distance:
                return {
                    'type': 'inverted',
                    'junction_pos': (primary['start'] + supp['start']) // 2,
                    'distance': abs(primary['start'] - supp['start']),
                    'orientation': 'inverted'
                }
        
        return None
    
    def _calculate_split_quality(self, primary, supp):
        """Calculate quality score for split alignment"""
        # Combine multiple quality metrics
        mapq_score = (primary['mapq'] + supp['mapq']) / 2 / 60  # Normalize to 0-1
        
        # Length score (longer alignments are more reliable)
        length_score = min(1.0, (primary['length'] + supp['length']) / 1000)
        
        # Edit distance score
        primary_nm_rate = primary['nm'] / primary['length'] if primary['length'] > 0 else 0
        supp_nm_rate = supp['nm'] / supp['length'] if supp['length'] > 0 else 0
        nm_score = 1 - (primary_nm_rate + supp_nm_rate) / 2
        
        # Combined score
        quality_score = (mapq_score + length_score + nm_score) / 3
        
        return quality_score
    
    def _create_split_candidates(self, validated_splits, chromosome):
        """Create candidates from validated split alignments"""
        # Group splits by junction position
        junction_groups = defaultdict(list)
        
        for split in validated_splits:
            junction_pos = split['pattern']['junction_pos']
            # Use broader binning for ONT data
            bin_pos = (junction_pos // 500) * 500
            junction_groups[bin_pos].append(split)
        
        candidates = []
        
        for bin_pos, splits in junction_groups.items():
            if len(splits) < self.min_support:
                continue
            
            # Calculate consensus junction position and boundaries
            junction_positions = [s['pattern']['junction_pos'] for s in splits]
            mean_junction = int(np.mean(junction_positions))
            
            # Estimate boundaries
            all_positions = []
            for split in splits:
                all_positions.extend([
                    split['primary']['start'],
                    split['primary']['end'],
                    split['supplementary']['start'],
                    split['supplementary']['end']
                ])
            
            if not all_positions:
                continue
            
            # Define circular DNA boundaries
            start = min(all_positions) - 200
            end = max(all_positions) + 200
            start = max(0, start)
            
            length = end - start
            
            # Filter by size
            if 200 <= length <= 100000:
                # Calculate support metrics
                unique_reads = len(set(s['read_name'] for s in splits))
                avg_quality = np.mean([s['quality_score'] for s in splits])
                
                # Pattern distribution
                pattern_types = defaultdict(int)
                for split in splits:
                    pattern_types[split['pattern']['type']] += 1
                
                # Create candidate
                candidate = CircularCandidate(
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    length=length,
                    split_support=unique_reads,
                    confidence_score=0.0,
                    detection_method='split_read'
                )
                
                # Add additional metrics
                candidate.split_quality = avg_quality
                candidate.split_patterns = dict(pattern_types)
                candidate.junction_position = mean_junction
                
                # Calculate confidence score
                candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
                
                # Apply quality-based confidence adjustment
                candidate.confidence_score *= avg_quality
                
                candidates.append(candidate)
        
        return candidates