#!/usr/bin/env python3
"""
CircONTrack module for characterizing circular DNA regions to determine their composition:
- Pure host (mouse) circles
- Pure viral circles  
- Chimeric host-viral constructs
- Integration sites
"""

import pysam
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
try:
    import pandas as pd
except ImportError:
    pd = None

# Import CircularCandidate from utils if available, otherwise define minimal version
try:
    from .utils import CircularCandidate
except ImportError:
    from dataclasses import dataclass
    
    @dataclass
    class CircularCandidate:
        """Minimal CircularCandidate for standalone use"""
        chromosome: str
        start: int
        end: int
        length: int
        confidence_score: float = 0.0
        detection_method: str = 'unknown'
        mean_coverage: Optional[float] = None
        fold_enrichment: Optional[float] = None
        coverage_uniformity: Optional[float] = None

class CircularRegionCharacterization:
    """Store characterization results for a circular DNA region"""
    def __init__(self, region_id: str, chromosome: str, start: int, end: int, length: int,
                 total_reads: int = 0, host_only_reads: int = 0, viral_only_reads: int = 0,
                 chimeric_reads: int = 0, unmapped_reads: int = 0,
                 host_coverage_percent: float = 0.0, viral_coverage_percent: float = 0.0,
                 has_integration_junction: bool = False, integration_breakpoints: List[Tuple[int, int]] = None,
                 junction_read_count: int = 0, region_type: str = 'unknown', confidence: float = 0.0,
                 gc_content: float = 0.0, mean_read_length: float = 0.0, viral_species: List[str] = None):
        self.region_id = region_id
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.length = length
        self.total_reads = total_reads
        self.host_only_reads = host_only_reads
        self.viral_only_reads = viral_only_reads
        self.chimeric_reads = chimeric_reads
        self.unmapped_reads = unmapped_reads
        self.host_coverage_percent = host_coverage_percent
        self.viral_coverage_percent = viral_coverage_percent
        self.has_integration_junction = has_integration_junction
        self.integration_breakpoints = integration_breakpoints or []
        self.junction_read_count = junction_read_count
        self.region_type = region_type
        self.confidence = confidence
        self.gc_content = gc_content
        self.mean_read_length = mean_read_length
        self.viral_species = viral_species or []

class CircularDNACharacterizer:
    """Characterize circular DNA regions for their composition"""
    
    def __init__(self, host_ref_path: str, combined_ref_path: str, 
                 viral_contig_prefix: str = "viral"):
        """
        Initialize with reference paths
        
        Args:
            host_ref_path: Path to host-only reference (mm10)
            combined_ref_path: Path to combined host+viral reference
            viral_contig_prefix: Prefix identifying viral contigs in combined reference
        """
        self.host_ref = pysam.FastaFile(host_ref_path)
        self.combined_ref = pysam.FastaFile(combined_ref_path)
        self.viral_prefix = viral_contig_prefix
        
        # Build viral contig list from combined reference
        self.viral_contigs = [ref for ref in self.combined_ref.references 
                              if ref.startswith(viral_contig_prefix)]
        self.host_contigs = [ref for ref in self.combined_ref.references 
                            if not ref.startswith(viral_contig_prefix)]
    
    def characterize_candidates(self, candidates: List[CircularCandidate],
                               host_bam_path: str, combined_bam_path: str,
                               min_mapq: int = 20) -> List[CircularRegionCharacterization]:
        """
        Characterize a list of CircularCandidate objects from CircONTrack
        
        Args:
            candidates: List of CircularCandidate objects from detection
            host_bam_path: BAM aligned to host-only reference
            combined_bam_path: BAM aligned to combined reference
            min_mapq: Minimum mapping quality threshold
            
        Returns:
            List of characterization results
        """
        characterizations = []
        
        for candidate in candidates:
            # Create BED-like string for compatibility
            bed_line = f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
            bed_line += f"circDNA_{len(characterizations)+1}\t{candidate.confidence_score}"
            
            char = self.characterize_region(bed_line, host_bam_path, combined_bam_path, min_mapq)
            
            # Add original detection method info
            char.detection_method = candidate.detection_method
            char.original_confidence = candidate.confidence_score
            
            characterizations.append(char)
        
        return characterizations
    
    def characterize_region(self, 
                           region_bed_line: str,
                           host_bam_path: str,
                           combined_bam_path: str,
                           min_mapq: int = 20) -> CircularRegionCharacterization:
        """
        Characterize a single circular DNA region
        
        Args:
            region_bed_line: BED format line with detected region
            host_bam_path: BAM aligned to host-only reference
            combined_bam_path: BAM aligned to combined reference
            min_mapq: Minimum mapping quality threshold
        """
        # Parse BED line
        fields = region_bed_line.strip().split('\t')
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        region_id = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
        length = end - start
        
        # Open BAM files
        host_bam = pysam.AlignmentFile(host_bam_path, "rb")
        combined_bam = pysam.AlignmentFile(combined_bam_path, "rb")
        
        # Analyze read composition
        read_analysis = self._analyze_read_composition(
            chrom, start, end, host_bam, combined_bam, min_mapq
        )
        
        # Analyze coverage patterns
        coverage_analysis = self._analyze_coverage_patterns(
            chrom, start, end, combined_bam
        )
        
        # Detect integration junctions
        integration_analysis = self._detect_integration_junctions(
            chrom, start, end, combined_bam, min_mapq
        )
        
        # Classify region type
        region_type, confidence = self._classify_region(
            read_analysis, coverage_analysis, integration_analysis
        )
        
        # Calculate additional metrics
        gc_content = self._calculate_gc_content(chrom, start, end)
        mean_read_length = read_analysis.get('mean_read_length', 0)
        viral_species = self._identify_viral_species(read_analysis)
        
        # Close BAM files
        host_bam.close()
        combined_bam.close()
        
        return CircularRegionCharacterization(
            region_id=region_id,
            chromosome=chrom,
            start=start,
            end=end,
            length=length,
            total_reads=read_analysis['total_reads'],
            host_only_reads=read_analysis['host_only'],
            viral_only_reads=read_analysis['viral_only'],
            chimeric_reads=read_analysis['chimeric'],
            unmapped_reads=read_analysis['unmapped'],
            host_coverage_percent=coverage_analysis['host_coverage_percent'],
            viral_coverage_percent=coverage_analysis['viral_coverage_percent'],
            has_integration_junction=integration_analysis['has_junction'],
            integration_breakpoints=integration_analysis['breakpoints'],
            junction_read_count=integration_analysis['junction_reads'],
            region_type=region_type,
            confidence=confidence,
            gc_content=gc_content,
            mean_read_length=mean_read_length,
            viral_species=viral_species
        )
    
    def _analyze_read_composition(self, chrom: str, start: int, end: int,
                                 host_bam: pysam.AlignmentFile,
                                 combined_bam: pysam.AlignmentFile,
                                 min_mapq: int) -> Dict:
        """Analyze the composition of reads in the region"""
        read_stats = {
            'total_reads': 0,
            'host_only': 0,
            'viral_only': 0,
            'chimeric': 0,
            'unmapped': 0,
            'read_lengths': [],
            'viral_refs': defaultdict(int)
        }
        
        # Get all reads from the region in combined BAM
        read_ids = set()
        for read in combined_bam.fetch(chrom, start, end):
            if read.mapping_quality >= min_mapq and not read.is_secondary:
                read_ids.add(read.query_name)
        
        # Analyze each unique read
        for read_id in read_ids:
            read_stats['total_reads'] += 1
            
            # Check if read maps to host in host-only BAM
            maps_to_host = self._read_maps_to_host(read_id, host_bam, min_mapq)
            
            # Check read's alignments in combined BAM
            alignments = self._get_read_alignments(read_id, combined_bam, min_mapq)
            
            # Classify read based on alignments
            has_host = any(not ref.startswith(self.viral_prefix) for ref in alignments['refs'])
            has_viral = any(ref.startswith(self.viral_prefix) for ref in alignments['refs'])
            
            if has_host and has_viral:
                read_stats['chimeric'] += 1
            elif has_viral and not maps_to_host:
                read_stats['viral_only'] += 1
                for ref in alignments['refs']:
                    if ref.startswith(self.viral_prefix):
                        read_stats['viral_refs'][ref] += 1
            elif has_host or maps_to_host:
                read_stats['host_only'] += 1
            else:
                read_stats['unmapped'] += 1
            
            if alignments['read_length']:
                read_stats['read_lengths'].append(alignments['read_length'])
        
        read_stats['mean_read_length'] = np.mean(read_stats['read_lengths']) if read_stats['read_lengths'] else 0
        
        return read_stats
    
    def _analyze_coverage_patterns(self, chrom: str, start: int, end: int,
                                   combined_bam: pysam.AlignmentFile) -> Dict:
        """Analyze coverage patterns to identify viral vs host regions"""
        window_size = 100
        num_windows = (end - start) // window_size + 1
        
        host_coverage = np.zeros(num_windows)
        viral_coverage = np.zeros(num_windows)
        
        for i, pos in enumerate(range(start, end, window_size)):
            window_end = min(pos + window_size, end)
            
            for pileup_col in combined_bam.pileup(chrom, pos, window_end, truncate=True):
                window_idx = (pileup_col.pos - start) // window_size
                if window_idx >= num_windows:
                    continue
                
                for pileup_read in pileup_col.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        read = pileup_read.alignment
                        
                        # Check if read has supplementary alignment to viral
                        if self._has_viral_alignment(read, combined_bam):
                            viral_coverage[window_idx] += 1
                        else:
                            host_coverage[window_idx] += 1
        
        # Calculate percentages
        total_coverage = host_coverage + viral_coverage
        host_percent = np.sum(host_coverage > 0) / num_windows * 100 if num_windows > 0 else 0
        viral_percent = np.sum(viral_coverage > 0) / num_windows * 100 if num_windows > 0 else 0
        
        return {
            'host_coverage_percent': host_percent,
            'viral_coverage_percent': viral_percent,
            'host_coverage_array': host_coverage,
            'viral_coverage_array': viral_coverage
        }
    
    def _detect_integration_junctions(self, chrom: str, start: int, end: int,
                                     combined_bam: pysam.AlignmentFile,
                                     min_mapq: int) -> Dict:
        """Detect potential integration junctions"""
        junctions = []
        junction_reads = 0
        
        # Look for split/chimeric reads that span host-viral boundary
        for read in combined_bam.fetch(chrom, start, end):
            if read.mapping_quality < min_mapq or read.is_secondary:
                continue
            
            # Check for SA tag (supplementary alignments)
            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                
                # Parse supplementary alignments
                for sa in sa_tag.rstrip(';').split(';'):
                    sa_fields = sa.split(',')
                    sa_ref = sa_fields[0]
                    sa_pos = int(sa_fields[1])
                    
                    # Check if supplementary alignment is to viral
                    if sa_ref.startswith(self.viral_prefix):
                        # This is a potential integration junction
                        junction_pos = read.reference_start if read.reference_start > start else read.reference_end
                        junctions.append((junction_pos, sa_pos))
                        junction_reads += 1
            
            # Also check for soft-clipped reads at region boundaries
            if read.cigartuples:
                if read.cigartuples[0][0] == 4:  # Soft clip at start
                    clip_len = read.cigartuples[0][1]
                    if clip_len > 20:  # Significant soft clip
                        junction_reads += 1
                        junctions.append((read.reference_start, -1))
                
                if read.cigartuples[-1][0] == 4:  # Soft clip at end
                    clip_len = read.cigartuples[-1][1]
                    if clip_len > 20:
                        junction_reads += 1
                        junctions.append((read.reference_end, -1))
        
        return {
            'has_junction': len(junctions) > 0,
            'breakpoints': junctions,
            'junction_reads': junction_reads
        }
    
    def _classify_region(self, read_analysis: Dict, coverage_analysis: Dict,
                        integration_analysis: Dict) -> Tuple[str, float]:
        """Classify the region based on all analyses"""
        total_reads = read_analysis['total_reads']
        if total_reads == 0:
            return 'unknown', 0.0
        
        # Calculate proportions
        host_prop = read_analysis['host_only'] / total_reads
        viral_prop = read_analysis['viral_only'] / total_reads
        chimeric_prop = read_analysis['chimeric'] / total_reads
        
        # Classification logic
        if integration_analysis['has_junction'] and chimeric_prop > 0.1:
            return 'integration_site', min(0.95, chimeric_prop + 0.3)
        
        elif viral_prop > 0.8:
            return 'viral_only', viral_prop
        
        elif host_prop > 0.8:
            return 'host_only', host_prop
        
        elif chimeric_prop > 0.3 or (viral_prop > 0.2 and host_prop > 0.2):
            return 'chimeric', max(chimeric_prop, (viral_prop + host_prop) / 2)
        
        else:
            return 'mixed', 0.5
    
    def _read_maps_to_host(self, read_id: str, host_bam: pysam.AlignmentFile,
                           min_mapq: int) -> bool:
        """Check if a read maps to host-only reference"""
        try:
            for read in host_bam.fetch(until_eof=True):
                if read.query_name == read_id and read.mapping_quality >= min_mapq:
                    return True
        except:
            pass
        return False
    
    def _get_read_alignments(self, read_id: str, bam: pysam.AlignmentFile,
                            min_mapq: int) -> Dict:
        """Get all alignments for a read"""
        alignments = {
            'refs': [],
            'positions': [],
            'read_length': None
        }
        
        for read in bam.fetch(until_eof=True):
            if read.query_name == read_id and read.mapping_quality >= min_mapq:
                alignments['refs'].append(read.reference_name)
                alignments['positions'].append(read.reference_start)
                if alignments['read_length'] is None:
                    alignments['read_length'] = read.query_length
        
        return alignments
    
    def _has_viral_alignment(self, read: pysam.AlignedSegment,
                            bam: pysam.AlignmentFile) -> bool:
        """Check if read has any alignment to viral sequences"""
        # Check SA tag for supplementary alignments
        if read.has_tag('SA'):
            sa_tag = read.get_tag('SA')
            for sa in sa_tag.rstrip(';').split(';'):
                sa_ref = sa.split(',')[0]
                if sa_ref.startswith(self.viral_prefix):
                    return True
        
        # Check if primary alignment is to viral
        if read.reference_name and read.reference_name.startswith(self.viral_prefix):
            return True
        
        return False
    
    def _calculate_gc_content(self, chrom: str, start: int, end: int) -> float:
        """Calculate GC content of the region"""
        try:
            # Try to get sequence from combined reference first
            if chrom in self.combined_ref.references:
                seq = self.combined_ref.fetch(chrom, start, end).upper()
            else:
                seq = self.host_ref.fetch(chrom, start, end).upper()
            
            gc_count = seq.count('G') + seq.count('C')
            return gc_count / len(seq) if len(seq) > 0 else 0.0
        except:
            return 0.0
    
    def _identify_viral_species(self, read_analysis: Dict) -> List[str]:
        """Identify viral species present in the region"""
        viral_species = []
        for viral_ref, count in read_analysis['viral_refs'].items():
            if count > 0:
                # Extract species name from reference name
                species = viral_ref.replace(self.viral_prefix, '').strip('_')
                viral_species.append(species)
        
        return viral_species
    
    def generate_summary_report(self, characterizations: List[CircularRegionCharacterization]):
        """Generate summary report for all characterized regions"""
        if pd is None:
            # Return dictionary if pandas not available
            return self._generate_dict_report(characterizations)
        
        data = []
        for char in characterizations:
            data.append({
                'Region_ID': char.region_id,
                'Chromosome': char.chromosome,
                'Start': char.start,
                'End': char.end,
                'Length': char.length,
                'Type': char.region_type,
                'Detection_Method': getattr(char, 'detection_method', 'unknown'),
                'Original_Confidence': getattr(char, 'original_confidence', 0.0),
                'Characterization_Confidence': f"{char.confidence:.3f}",
                'Total_Reads': char.total_reads,
                'Host_Only_Reads': char.host_only_reads,
                'Viral_Only_Reads': char.viral_only_reads,
                'Chimeric_Reads': char.chimeric_reads,
                'Host_Coverage_%': f"{char.host_coverage_percent:.1f}",
                'Viral_Coverage_%': f"{char.viral_coverage_percent:.1f}",
                'Has_Integration': char.has_integration_junction,
                'Junction_Reads': char.junction_read_count,
                'GC_Content': f"{char.gc_content:.3f}",
                'Mean_Read_Length': f"{char.mean_read_length:.0f}",
                'Viral_Species': ','.join(char.viral_species) if char.viral_species else 'None'
            })
        
        return pd.DataFrame(data)
    
    def _generate_dict_report(self, characterizations: List[CircularRegionCharacterization]) -> List[Dict]:
        """Generate report as list of dictionaries when pandas not available"""
        data = []
        for char in characterizations:
            data.append({
                'Region_ID': char.region_id,
                'Chromosome': char.chromosome,
                'Start': char.start,
                'End': char.end,
                'Length': char.length,
                'Type': char.region_type,
                'Detection_Method': getattr(char, 'detection_method', 'unknown'),
                'Original_Confidence': getattr(char, 'original_confidence', 0.0),
                'Characterization_Confidence': char.confidence,
                'Total_Reads': char.total_reads,
                'Host_Only_Reads': char.host_only_reads,
                'Viral_Only_Reads': char.viral_only_reads,
                'Chimeric_Reads': char.chimeric_reads,
                'Host_Coverage_%': char.host_coverage_percent,
                'Viral_Coverage_%': char.viral_coverage_percent,
                'Has_Integration': char.has_integration_junction,
                'Junction_Reads': char.junction_read_count,
                'GC_Content': char.gc_content,
                'Mean_Read_Length': char.mean_read_length,
                'Viral_Species': ','.join(char.viral_species) if char.viral_species else 'None'
            })
        return data
    
    def write_characterized_bed(self, characterizations: List[CircularRegionCharacterization],
                               output_file: str):
        """Write characterized regions to enhanced BED format"""
        with open(output_file, 'w') as f:
            # Write header
            f.write("#chr\tstart\tend\tname\tscore\tstrand\tdetection_method\tlength\t")
            f.write("region_type\tchar_confidence\thost_reads\tviral_reads\tchimeric_reads\t")
            f.write("has_integration\tviral_species\n")
            
            for char in characterizations:
                # Convert confidence to BED score (0-1000)
                bed_score = int(char.confidence * 1000)
                
                # Format viral species
                viral_str = ','.join(char.viral_species) if char.viral_species else 'none'
                
                # Write BED line with extended fields
                f.write(f"{char.chromosome}\t{char.start}\t{char.end}\t")
                f.write(f"{char.region_id}\t{bed_score}\t.\t")
                f.write(f"{getattr(char, 'detection_method', 'unknown')}\t{char.length}\t")
                f.write(f"{char.region_type}\t{char.confidence:.3f}\t")
                f.write(f"{char.host_only_reads}\t{char.viral_only_reads}\t{char.chimeric_reads}\t")
                f.write(f"{char.has_integration_junction}\t{viral_str}\n")


def integrate_with_circontrack(detector, bam_file: str, reference_fasta: str,
                              host_ref: str, combined_ref: str, 
                              host_bam: str, combined_bam: str,
                              viral_prefix: str = "viral",
                              output_prefix: str = "circontrack"):
    """
    Integration function to use characterization with CircONTrack pipeline
    
    Args:
        detector: CircularDNADetector instance from CircONTrack
        bam_file: Original BAM file for CircONTrack detection
        reference_fasta: Reference FASTA for CircONTrack
        host_ref: Host-only reference FASTA
        combined_ref: Combined host+viral reference FASTA
        host_bam: BAM aligned to host-only reference
        combined_bam: BAM aligned to combined reference
        viral_prefix: Prefix for viral contigs
        output_prefix: Prefix for output files
    
    Returns:
        Tuple of (detected_candidates, characterizations)
    """
    print("Step 1: Detecting circular DNA with CircONTrack...")
    candidates = detector.detect_circular_dna(bam_file, reference_fasta)
    print(f"  Found {len(candidates)} circular DNA candidates")
    
    print("\nStep 2: Characterizing regions for viral content...")
    characterizer = CircularDNACharacterizer(host_ref, combined_ref, viral_prefix)
    characterizations = characterizer.characterize_candidates(
        candidates, host_bam, combined_bam
    )
    
    # Update candidates with characterization info
    for candidate, char in zip(candidates, characterizations):
        candidate.region_type = char.region_type
        candidate.has_viral = char.viral_only_reads > 0 or char.chimeric_reads > 0
        candidate.is_integration = char.has_integration_junction
        candidate.viral_species = char.viral_species
    
    # Write outputs
    print(f"\nStep 3: Writing results to {output_prefix}_characterized.bed")
    characterizer.write_characterized_bed(characterizations, f"{output_prefix}_characterized.bed")
    
    # Generate summary report
    if pd is not None:
        summary_df = characterizer.generate_summary_report(characterizations)
        summary_df.to_csv(f"{output_prefix}_summary.tsv", sep='\t', index=False)
        print(f"  Summary report saved to {output_prefix}_summary.tsv")
        
        # Print statistics
        print("\n=== Detection Summary ===")
        print(f"Total circular DNA regions: {len(candidates)}")
        type_counts = summary_df['Type'].value_counts()
        for region_type, count in type_counts.items():
            print(f"  {region_type}: {count} ({count/len(candidates)*100:.1f}%)")
    
    return candidates, characterizations


def main():
    """Main execution function for standalone or integrated use"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Characterize circular DNA regions for viral content')
    
    # Add subcommands for different modes
    subparsers = parser.add_subparsers(dest='mode', help='Operation mode')
    
    # Standalone mode - characterize existing BED file
    standalone = subparsers.add_parser('characterize', help='Characterize regions from BED file')
    standalone.add_argument('bed_file', help='BED file with CircONTrack results')
    standalone.add_argument('host_ref', help='Host-only reference FASTA (mm10)')
    standalone.add_argument('combined_ref', help='Combined host+viral reference FASTA')
    standalone.add_argument('host_bam', help='BAM aligned to host-only reference')
    standalone.add_argument('combined_bam', help='BAM aligned to combined reference')
    standalone.add_argument('-o', '--output', default='circdna_characterization',
                          help='Output prefix')
    standalone.add_argument('--viral-prefix', default='viral',
                          help='Prefix for viral contigs in combined reference')
    standalone.add_argument('--min-mapq', type=int, default=20,
                          help='Minimum mapping quality')
    
    # Integrated mode - run full CircONTrack + characterization
    integrated = subparsers.add_parser('full', help='Run full CircONTrack detection + characterization')
    integrated.add_argument('bam_file', help='Input BAM file for CircONTrack')
    integrated.add_argument('reference', help='Reference FASTA for CircONTrack')
    integrated.add_argument('host_ref', help='Host-only reference FASTA')
    integrated.add_argument('combined_ref', help='Combined host+viral reference FASTA')
    integrated.add_argument('host_bam', help='BAM aligned to host-only reference')
    integrated.add_argument('combined_bam', help='BAM aligned to combined reference')
    integrated.add_argument('-o', '--output', default='circontrack',
                          help='Output prefix')
    integrated.add_argument('--viral-prefix', default='viral',
                          help='Prefix for viral contigs')
    integrated.add_argument('--min-fold-enrichment', type=float, default=1.5,
                          help='Minimum fold enrichment for CircONTrack')
    integrated.add_argument('--min-coverage', type=int, default=5,
                          help='Minimum coverage for CircONTrack')
    
    args = parser.parse_args()
    
    if not args.mode:
        parser.print_help()
        return
    
    if args.mode == 'characterize':
        # Standalone characterization mode
        characterizer = CircularDNACharacterizer(
            args.host_ref,
            args.combined_ref,
            args.viral_prefix
        )
        
        # Process each region from BED file
        characterizations = []
        with open(args.bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                print(f"Processing: {line.strip()}")
                char = characterizer.characterize_region(
                    line,
                    args.host_bam,
                    args.combined_bam,
                    args.min_mapq
                )
                characterizations.append(char)
                
                # Print immediate summary
                print(f"  -> Type: {char.region_type} (confidence: {char.confidence:.3f})")
                print(f"     Reads: {char.host_only_reads} host, {char.viral_only_reads} viral, "
                      f"{char.chimeric_reads} chimeric")
                if char.viral_species:
                    print(f"     Viral species: {', '.join(char.viral_species)}")
                if char.has_integration_junction:
                    print(f"     Integration detected with {char.junction_read_count} junction reads")
                print()
        
        # Write outputs
        characterizer.write_characterized_bed(characterizations, f"{args.output}_characterized.bed")
        
        if pd is not None:
            summary_df = characterizer.generate_summary_report(characterizations)
            summary_df.to_csv(f"{args.output}_summary.tsv", sep='\t', index=False)
            print(f"Summary report saved to: {args.output}_summary.tsv")
    
    elif args.mode == 'full':
        # Integrated mode with CircONTrack
        try:
            from circDNA_detection import CircularDNADetector
            
            # Initialize CircONTrack detector
            detector = CircularDNADetector(
                min_fold_enrichment=args.min_fold_enrichment,
                min_coverage=args.min_coverage
            )
            
            # Run integrated pipeline
            candidates, characterizations = integrate_with_circontrack(
                detector,
                args.bam_file,
                args.reference,
                args.host_ref,
                args.combined_ref,
                args.host_bam,
                args.combined_bam,
                args.viral_prefix,
                args.output
            )
            
        except ImportError:
            print("Error: CircONTrack not found. Please install CircONTrack first.")
            print("For standalone characterization, use 'characterize' mode instead.")
            return


if __name__ == "__main__":
    main()