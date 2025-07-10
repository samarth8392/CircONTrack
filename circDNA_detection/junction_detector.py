#!/usr/bin/env python3
"""
Circular DNA Junction Detection
Detects back-to-back junctions characteristic of circular DNA
"""

import pysam
import numpy as np
from collections import defaultdict, Counter
from .utils import CircularCandidate

class JunctionDetector:
    def __init__(self, min_support=3, max_junction_distance=1000):
        self.min_support = min_support
        self.max_junction_distance = max_junction_distance
    
    def detect_junctions(self, bam_file, chromosome=None):
        """Detect circular DNA junctions from discordant read pairs"""
        print("  Detecting junctions...")
        
        junction_candidates = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chrom in chromosomes:
                junctions = self._find_chromosome_junctions(bam, chrom)
                candidates = self._cluster_junctions(junctions, chrom)
                junction_candidates.extend(candidates)
        
        return junction_candidates
    
    def _find_chromosome_junctions(self, bam, chromosome):
        """Find potential circular junctions in chromosome"""
        junctions = defaultdict(list)
        
        for read in bam.fetch(chromosome):
            if (read.is_unmapped or read.mate_is_unmapped or 
                read.is_secondary or read.is_supplementary):
                continue
            
            # Look for discordant pairs (back-to-back orientation)
            if self._is_circular_junction_pair(read):
                # Record junction point
                if read.is_reverse:
                    junction_pos = read.reference_start
                else:
                    junction_pos = read.reference_end
                
                junctions[junction_pos].append(read)
        
        return junctions
    
    def _is_circular_junction_pair(self, read):
        """Check if read pair indicates circular DNA junction"""
        # Require proper pair mapping to same chromosome
        if (read.reference_name != read.next_reference_name or
            not read.is_proper_pair):
            return False
        
        # Check for back-to-back orientation pattern
        # Forward read followed by reverse read (or vice versa)
        # with small insert size indicating circular junction
        
        insert_size = abs(read.template_length)
        
        # Back-to-back reads with small insert size
        if insert_size < self.max_junction_distance:
            # Check orientation pattern
            if read.is_reverse != read.mate_is_reverse:
                return True
        
        return False
    
    def _cluster_junctions(self, junctions, chromosome):
        """Cluster junction points to identify circular DNA regions"""
        candidates = []
        
        if not junctions:
            return candidates
        
        # Sort junction positions
        junction_positions = sorted(junctions.keys())
        
        # Cluster nearby junctions
        clusters = self._cluster_positions(junction_positions, max_distance=500)
        
        for cluster in clusters:
            if len(cluster) < 2:  # Need at least 2 junction points
                continue
            
            # Calculate cluster statistics
            cluster_start = min(cluster)
            cluster_end = max(cluster)
            
            # Count supporting reads
            supporting_reads = []
            for pos in cluster:
                supporting_reads.extend(junctions[pos])
            
            junction_support = len(supporting_reads)
            
            if junction_support < self.min_support:
                continue
            
            # Estimate circular DNA boundaries
            # Look for read pairs that span the potential circular region
            boundaries = self._estimate_circular_boundaries(
                supporting_reads, cluster_start, cluster_end
            )
            
            if boundaries:
                start, end = boundaries
                length = end - start
                
                # Filter by reasonable size
                if 200 <= length <= 100000:
                    candidate = CircularCandidate(
                        chromosome=chromosome,
                        start=start,
                        end=end,
                        length=length,
                        junction_support=junction_support,
                        confidence_score=0.0,
                        detection_method='junction'
                    )
                    candidates.append(candidate)
        
        return candidates
    
    def _cluster_positions(self, positions, max_distance=500):
        """Cluster nearby positions"""
        if not positions:
            return []
        
        clusters = []
        current_cluster = [positions[0]]
        
        for pos in positions[1:]:
            if pos - current_cluster[-1] <= max_distance:
                current_cluster.append(pos)
            else:
                clusters.append(current_cluster)
                current_cluster = [pos]
        
        clusters.append(current_cluster)
        return clusters
    
    def _estimate_circular_boundaries(self, reads, junction_start, junction_end):
        """Estimate circular DNA boundaries from supporting reads"""
        # Collect read positions
        read_starts = []
        read_ends = []
        
        for read in reads:
            read_starts.append(read.reference_start)
            read_ends.append(read.reference_end)
            
            # Also consider mate positions
            if read.next_reference_start:
                read_starts.append(read.next_reference_start)
                read_ends.append(read.next_reference_start + read.query_length)
        
        if not read_starts:
            return None
        
        # Estimate boundaries as outer bounds of supporting reads
        estimated_start = min(read_starts)
        estimated_end = max(read_ends)
        
        # Extend slightly to capture full circular element
        padding = 100
        estimated_start = max(0, estimated_start - padding)
        estimated_end = estimated_end + padding
        
        return estimated_start, estimated_end