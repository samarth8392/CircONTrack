#!/usr/bin/env python3
"""
Split Read Analysis for Circular DNA Detection
Analyzes supplementary alignments to identify circular DNA junctions
"""

import pysam
import numpy as np
from collections import defaultdict
from .utils import CircularCandidate


class SplitReadAnalyzer:
    """Analyzes split reads to detect circular DNA junctions"""
    
    def __init__(self, min_split_length=50, min_support=3, max_distance=1000):
        self.min_split_length = min_split_length
        self.min_support = min_support
        self.max_distance = max_distance
    
    def analyze_split_reads(self, bam_file, chromosome=None):
        """Analyze split reads to detect circular DNA junctions"""
        print("  Analyzing split reads...")
        
        split_candidates = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chrom in chromosomes:
                splits = self._collect_split_alignments(bam, chrom)
                candidates = self._analyze_split_patterns(splits, chrom)
                split_candidates.extend(candidates)
        
        return split_candidates
    
    def _collect_split_alignments(self, bam, chromosome):
        """Collect split read alignments for analysis"""
        split_reads = defaultdict(list)
        
        for read in bam.fetch(chromosome):
            if read.is_unmapped or read.is_secondary:
                continue
            
            # Look for supplementary alignments (split reads)
            if read.is_supplementary:
                continue
            
            # Check if read has supplementary alignments
            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                supplementary_alns = self._parse_sa_tag(sa_tag)
                
                # Filter supplementary alignments
                valid_supplements = []
                for supp in supplementary_alns:
                    if (supp['length'] >= self.min_split_length and
                        supp['chromosome'] == chromosome):
                        valid_supplements.append(supp)
                
                if valid_supplements:
                    split_reads[read.query_name].append({
                        'primary': {
                            'chromosome': chromosome,
                            'start': read.reference_start,
                            'end': read.reference_end,
                            'strand': '-' if read.is_reverse else '+',
                            'mapq': read.mapping_quality
                        },
                        'supplementary': valid_supplements
                    })
        
        return split_reads
    
    def _parse_sa_tag(self, sa_tag):
        """Parse SA tag to extract supplementary alignments"""
        supplements = []
        
        for alignment in sa_tag.split(';'):
            if not alignment:
                continue
            
            parts = alignment.split(',')
            if len(parts) >= 6:
                chromosome = parts[0]
                pos = int(parts[1])
                strand = parts[2]
                cigar = parts[3]
                mapq = int(parts[4])
                
                # Estimate alignment length from CIGAR
                length = self._estimate_alignment_length(cigar)
                
                supplements.append({
                    'chromosome': chromosome,
                    'position': pos,
                    'strand': strand,
                    'cigar': cigar,
                    'mapq': mapq,
                    'length': length
                })
        
        return supplements
    
    def _estimate_alignment_length(self, cigar_string):
        """Estimate alignment length from CIGAR string"""
        import re
        
        # Extract match/mismatch operations
        matches = re.findall(r'(\d+)[M=X]', cigar_string)
        return sum(int(match) for match in matches)
    
    def _analyze_split_patterns(self, split_reads, chromosome):
        """Analyze split read patterns to identify circular DNA"""
        candidates = []
        
        # Group split reads by junction points
        junction_groups = defaultdict(list)
        
        for read_name, split_info_list in split_reads.items():
            for split_info in split_info_list:
                primary = split_info['primary']
                
                for supp in split_info['supplementary']:
                    # Calculate potential junction points
                    if primary['strand'] == '+':
                        junction1 = primary['end']
                    else:
                        junction1 = primary['start']
                    
                    if supp['strand'] == '+':
                        junction2 = supp['position']
                    else:
                        junction2 = supp['position'] + supp['length']
                    
                    # Check for circular pattern (back-to-back alignment)
                    if abs(junction1 - junction2) <= self.max_distance:
                        junction_key = (min(junction1, junction2), max(junction1, junction2))
                        junction_groups[junction_key].append({
                            'read_name': read_name,
                            'primary': primary,
                            'supplementary': supp,
                            'junction1': junction1,
                            'junction2': junction2
                        })
        
        # Analyze junction groups
        for junction_key, supporting_reads in junction_groups.items():
            if len(supporting_reads) < self.min_support:
                continue
            
            # Estimate circular DNA boundaries
            boundaries = self._estimate_boundaries_from_splits(supporting_reads)
            
            if boundaries:
                start, end = boundaries
                length = end - start
                
                if 200 <= length <= 100000:
                    # Create candidate WITHOUT confidence score
                    candidate = CircularCandidate(
                        chromosome=chromosome,
                        start=start,
                        end=end,
                        length=length,
                        split_support=len(supporting_reads),
                        confidence_score=0.0,  # Placeholder - will be set by ConfidenceScorer
                        detection_method='split_read'
                    )
                    candidates.append(candidate)
        
        return candidates
    
    def _estimate_boundaries_from_splits(self, supporting_reads):
        """Estimate circular DNA boundaries from split read alignments"""
        all_positions = []
        
        for read_info in supporting_reads:
            primary = read_info['primary']
            supp = read_info['supplementary']
            
            all_positions.extend([
                primary['start'], primary['end'],
                supp['position'], supp['position'] + supp['length']
            ])
        
        if not all_positions:
            return None
        
        # Use outer bounds with small padding
        start = min(all_positions) - 50
        end = max(all_positions) + 50
        
        return max(0, start), end