#!/usr/bin/env python3
"""
ONT-Optimized Circular DNA Junction Detection
Detects circular DNA junctions from long-read sequencing data
"""

import pysam
import numpy as np
from collections import defaultdict
import re
from .utils import CircularCandidate
from .confidence_scorer import ConfidenceScorer

class JunctionDetector:
    def __init__(self, min_support=3, max_junction_distance=1000, 
             min_clip_length=50, min_alignment_length=100,
             min_mapping_quality=20, verbose=False):
        """
        Initialize junction detector for ONT data
        
        Parameters:
        -----------
        min_support : int
            Minimum number of reads supporting a junction
        max_junction_distance : int
            Maximum distance between junction endpoints to consider them the same
        min_clip_length : int
            Minimum soft-clip length to consider for junction detection
        min_alignment_length : int
            Minimum alignment length for supplementary alignments
        min_mapping_quality : int
            Minimum mapping quality for reads
        """
        self.min_support = min_support
        self.max_junction_distance = max_junction_distance
        self.min_clip_length = min_clip_length
        self.min_alignment_length = min_alignment_length
        self.min_mapping_quality = min_mapping_quality
        self.verbose = verbose
        self.confidence_scorer = ConfidenceScorer()
    
    def detect_junctions(self, bam_file, chromosome=None):
        """Detect circular DNA junctions from ONT reads"""
        print("  Detecting junctions from long reads...")
        
        junction_candidates = []
        if self.verbose:
            print("  Detecting junctions from long reads...")
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            chromosomes = [chromosome] if chromosome else bam.references
            
            for chrom in chromosomes:
                # Collect junctions from different detection methods
                sa_junctions = self._find_sa_junctions(bam, chrom)
                clip_junctions = self._find_softclip_junctions(bam, chrom)
                chimeric_junctions = self._find_chimeric_junctions(bam, chrom)
                
                # Merge all junction evidence
                all_junctions = self._merge_junction_evidence(
                    sa_junctions, clip_junctions, chimeric_junctions, chrom
                )
                
                # Create candidates from clustered junctions
                candidates = self._create_junction_candidates(all_junctions, chrom)
                junction_candidates.extend(candidates)
        
        return junction_candidates
    
    def _find_sa_junctions(self, bam, chromosome):
        """Find junctions from supplementary alignments (SA tag)"""
        junctions = defaultdict(list)
        
        for read in bam.fetch(chromosome):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            if read.mapping_quality < self.min_mapping_quality:
                continue
            
            # Check for SA tag
            if not read.has_tag('SA'):
                continue
            
            sa_tag = read.get_tag('SA')
            sa_alignments = self._parse_sa_tag(sa_tag)
            
            # Analyze SA alignments for circular patterns
            primary_info = {
                'chr': chromosome,
                'start': read.reference_start,
                'end': read.reference_end,
                'strand': '-' if read.is_reverse else '+',
                'mapq': read.mapping_quality
            }
            
            for sa in sa_alignments:
                if sa['chr'] != chromosome:
                    continue
                
                if sa['mapq'] < self.min_mapping_quality:
                    continue
                
                # Check for circular junction pattern
                junction = self._check_circular_pattern(primary_info, sa)
                if junction:
                    junction['read_name'] = read.query_name
                    junction['type'] = 'SA'
                    junctions[junction['position']].append(junction)
        
        return junctions
    
    def _find_softclip_junctions(self, bam, chromosome):
        """Find junctions from soft-clipped reads"""
        junctions = defaultdict(list)
        
        for read in bam.fetch(chromosome):
            if read.is_unmapped or read.mapping_quality < self.min_mapping_quality:
                continue
            
            if not read.cigartuples:
                continue
            
            # Check for significant soft-clipping
            cigar = read.cigartuples
            
            # Left soft-clip
            left_clip = 0
            if cigar[0][0] == 4:  # 4 = soft clip
                left_clip = cigar[0][1]
            
            # Right soft-clip
            right_clip = 0
            if cigar[-1][0] == 4:
                right_clip = cigar[-1][1]
            
            # Check for junction patterns
            if left_clip >= self.min_clip_length:
                # Potential junction at read start
                junction = {
                    'position': read.reference_start,
                    'clip_length': left_clip,
                    'clip_side': 'left',
                    'read_name': read.query_name,
                    'type': 'soft_clip',
                    'strand': '-' if read.is_reverse else '+'
                }
                junctions[read.reference_start].append(junction)
            
            if right_clip >= self.min_clip_length:
                # Potential junction at read end
                junction = {
                    'position': read.reference_end,
                    'clip_length': right_clip,
                    'clip_side': 'right',
                    'read_name': read.query_name,
                    'type': 'soft_clip',
                    'strand': '-' if read.is_reverse else '+'
                }
                junctions[read.reference_end].append(junction)
            
            # Check for reads with large clips on both ends (junction-spanning reads)
            if left_clip >= self.min_clip_length and right_clip >= self.min_clip_length:
                # This read might span a junction
                junction = {
                    'position': (read.reference_start + read.reference_end) // 2,
                    'span_start': read.reference_start,
                    'span_end': read.reference_end,
                    'left_clip': left_clip,
                    'right_clip': right_clip,
                    'read_name': read.query_name,
                    'type': 'double_clip',
                    'strand': '-' if read.is_reverse else '+'
                }
                junctions[(read.reference_start + read.reference_end) // 2].append(junction)
        
        return junctions
    
    def _find_chimeric_junctions(self, bam, chromosome):
        """Find junctions from chimeric alignments (reads with multiple alignment blocks)"""
        junctions = defaultdict(list)
        
        for read in bam.fetch(chromosome):
            if read.is_unmapped or read.mapping_quality < self.min_mapping_quality:
                continue
            
            # Look for complex CIGAR strings that might indicate chimeric alignments
            if not read.cigartuples:
                continue
            
            # Check for patterns like: Match-LargeGap-Match
            cigar = read.cigartuples
            
            # Look for large deletions or reference skips that might indicate junctions
            current_pos = read.reference_start
            
            for op, length in cigar:
                if op == 0:  # Match
                    current_pos += length
                elif op == 2:  # Deletion
                    if length > self.max_junction_distance:
                        # Large deletion might indicate circular junction
                        junction = {
                            'position': current_pos,
                            'gap_size': length,
                            'read_name': read.query_name,
                            'type': 'chimeric_gap',
                            'strand': '-' if read.is_reverse else '+'
                        }
                        junctions[current_pos].append(junction)
                    current_pos += length
                elif op == 3:  # Reference skip
                    if length > self.max_junction_distance:
                        junction = {
                            'position': current_pos,
                            'skip_size': length,
                            'read_name': read.query_name,
                            'type': 'chimeric_skip',
                            'strand': '-' if read.is_reverse else '+'
                        }
                        junctions[current_pos].append(junction)
                    current_pos += length
        
        return junctions
    
    def _parse_sa_tag(self, sa_tag):
        """Parse supplementary alignment tag"""
        alignments = []
        
        for aln in sa_tag.strip().split(';'):
            if not aln:
                continue
            
            parts = aln.split(',')
            if len(parts) >= 6:
                chr_name = parts[0]
                position = int(parts[1]) - 1  # Convert to 0-based
                strand = parts[2]
                cigar = parts[3]
                mapq = int(parts[4])
                nm = int(parts[5]) if parts[5] else 0
                
                # Calculate alignment length from CIGAR
                aln_length = self._get_alignment_length(cigar)
                
                if aln_length >= self.min_alignment_length:
                    alignments.append({
                        'chr': chr_name,
                        'start': position,
                        'end': position + aln_length,
                        'strand': strand,
                        'cigar': cigar,
                        'mapq': mapq,
                        'nm': nm,
                        'length': aln_length
                    })
        
        return alignments
    
    def _get_alignment_length(self, cigar_string):
        """Calculate reference alignment length from CIGAR string"""
        # Parse CIGAR string
        length = 0
        matches = re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)
        
        for count, op in matches:
            count = int(count)
            if op in 'MDN=X':  # Operations that consume reference
                length += count
        
        return length
    
    def _check_circular_pattern(self, primary, supplementary):
        """Check if primary and supplementary alignments suggest circular junction"""
        # Look for back-to-back junction patterns
        
        # Case 1: Primary end connects to supplementary start
        if abs(primary['end'] - supplementary['start']) < self.max_junction_distance:
            return {
                'position': (primary['end'] + supplementary['start']) // 2,
                'junction_type': 'end_to_start',
                'support_type': 'SA',
                'primary_pos': primary['end'],
                'supp_pos': supplementary['start'],
                'distance': abs(primary['end'] - supplementary['start'])
            }
        
        # Case 2: Primary start connects to supplementary end
        if abs(primary['start'] - supplementary['end']) < self.max_junction_distance:
            return {
                'position': (primary['start'] + supplementary['end']) // 2,
                'junction_type': 'start_to_end',
                'support_type': 'SA',
                'primary_pos': primary['start'],
                'supp_pos': supplementary['end'],
                'distance': abs(primary['start'] - supplementary['end'])
            }
        
        # Case 3: Inverted orientation suggesting circular junction
        if primary['strand'] != supplementary['strand']:
            # Check for proximity of ends
            if abs(primary['end'] - supplementary['end']) < self.max_junction_distance:
                return {
                    'position': (primary['end'] + supplementary['end']) // 2,
                    'junction_type': 'inverted_ends',
                    'support_type': 'SA',
                    'primary_pos': primary['end'],
                    'supp_pos': supplementary['end'],
                    'distance': abs(primary['end'] - supplementary['end'])
                }
            
            if abs(primary['start'] - supplementary['start']) < self.max_junction_distance:
                return {
                    'position': (primary['start'] + supplementary['start']) // 2,
                    'junction_type': 'inverted_starts',
                    'support_type': 'SA',
                    'primary_pos': primary['start'],
                    'supp_pos': supplementary['start'],
                    'distance': abs(primary['start'] - supplementary['start'])
                }
        
        return None
    
    def _merge_junction_evidence(self, sa_junctions, clip_junctions, chimeric_junctions, chromosome):
        """Merge junction evidence from different sources"""
        all_junctions = defaultdict(list)
        
        # Combine all junction dictionaries
        for pos, junctions in sa_junctions.items():
            all_junctions[pos].extend(junctions)
        
        for pos, junctions in clip_junctions.items():
            all_junctions[pos].extend(junctions)
        
        for pos, junctions in chimeric_junctions.items():
            all_junctions[pos].extend(junctions)
        
        # Cluster nearby junctions
        clustered_junctions = self._cluster_junctions(all_junctions)
        
        return clustered_junctions
    
    def _cluster_junctions(self, junctions_dict):
        """Cluster nearby junction positions"""
        if not junctions_dict:
            return []
        
        # Sort positions
        positions = sorted(junctions_dict.keys())
        
        clusters = []
        current_cluster = {
            'positions': [positions[0]],
            'junctions': junctions_dict[positions[0]],
            'mean_pos': positions[0]
        }
        
        for pos in positions[1:]:
            if pos - current_cluster['mean_pos'] <= self.max_junction_distance:
                # Add to current cluster
                current_cluster['positions'].append(pos)
                current_cluster['junctions'].extend(junctions_dict[pos])
                current_cluster['mean_pos'] = int(np.mean(current_cluster['positions']))
            else:
                # Start new cluster
                clusters.append(current_cluster)
                current_cluster = {
                    'positions': [pos],
                    'junctions': junctions_dict[pos],
                    'mean_pos': pos
                }
        
        clusters.append(current_cluster)
        return clusters
    
    def _create_junction_candidates(self, junction_clusters, chromosome):
        """Create CircularCandidate objects from junction clusters"""
        candidates = []
        
        for cluster in junction_clusters:
            # Count support by type
            support_by_type = defaultdict(int)
            supporting_reads = set()
            
            for junction in cluster['junctions']:
                support_by_type[junction['type']] += 1
                supporting_reads.add(junction['read_name'])
            
            total_support = len(supporting_reads)
            
            if total_support < self.min_support:
                continue
            
            # Estimate circular DNA boundaries
            boundaries = self._estimate_boundaries_from_junctions(cluster, chromosome)
            
            if boundaries:
                start, end = boundaries
                length = end - start
                
                # Filter by reasonable size
                if 200 <= length <= 100000:
                    # Create candidate
                    candidate = CircularCandidate(
                        chromosome=chromosome,
                        start=start,
                        end=end,
                        length=length,
                        junction_support=total_support,
                        confidence_score=0.0,
                        detection_method='junction'
                    )
                    
                    # Add additional junction information
                    candidate.junction_types = dict(support_by_type)
                    candidate.junction_position = cluster['mean_pos']
                    
                    # Calculate confidence score
                    candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
                    
                    candidates.append(candidate)
        
        return candidates
    
    def _estimate_boundaries_from_junctions(self, junction_cluster, chromosome):
        """Estimate circular DNA boundaries from junction evidence"""
        # Collect all relevant positions
        positions = []
        
        for junction in junction_cluster['junctions']:
            positions.append(junction['position'])
            
            # Add span information if available
            if 'span_start' in junction:
                positions.append(junction['span_start'])
                positions.append(junction['span_end'])
            
            # Add primary/supplementary positions
            if 'primary_pos' in junction:
                positions.append(junction['primary_pos'])
            if 'supp_pos' in junction:
                positions.append(junction['supp_pos'])
        
        if not positions:
            return None
        
        # Use the range of positions with some padding
        min_pos = min(positions)
        max_pos = max(positions)
        
        # Extend to capture potential circular element
        padding = 500  # Larger padding for ONT
        start = max(0, min_pos - padding)
        end = max_pos + padding
        
        return start, end