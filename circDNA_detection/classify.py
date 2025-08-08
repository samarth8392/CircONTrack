#!/usr/bin/env python3
"""
CircONTrack Classify Module - Classify circular DNA candidates as host/viral/chimeric
Usage: circontrack classify <bed_file> <options>
"""

import pysam
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
import argparse
import sys
import os

# Import CircONTrack utilities if available
try:
    from utils import CircularCandidate
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

class CircularDNAClassifier:
    """Classify circular DNA regions by their genomic composition (host/viral/chimeric)"""
    
    def __init__(self, host_ref_path: str, combined_ref_path: str, 
                 viral_contig_prefix: str = "viral"):
        """
        Initialize classifier with reference genomes
        
        Args:
            host_ref_path: Path to host-only reference (e.g., mm10)
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
        
        print(f"Initialized classifier with {len(self.host_contigs)} host and {len(self.viral_contigs)} viral contigs")
    
    def classify_bed_file(self, bed_file: str, host_bam_path: str, 
                          combined_bam_path: str, output_prefix: str,
                          min_mapq: int = 20, verbose: bool = True):
        """
        Classify all regions from a BED file
        
        Args:
            bed_file: CircONTrack output BED file
            host_bam_path: BAM aligned to host-only reference
            combined_bam_path: BAM aligned to combined reference
            output_prefix: Prefix for output files
            min_mapq: Minimum mapping quality threshold
            verbose: Print progress messages
        """
        # Open BAM files
        host_bam = pysam.AlignmentFile(host_bam_path, "rb")
        combined_bam = pysam.AlignmentFile(combined_bam_path, "rb")
        
        classifications = []
        
        # Parse BED file
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3] if len(fields) > 3 else f"region_{line_num}"
                
                if verbose:
                    print(f"Classifying {name}: {chrom}:{start}-{end}")
                
                # Classify this region
                classification = self._classify_region(
                    chrom, start, end, name, host_bam, combined_bam, min_mapq
                )
                classifications.append(classification)
                
                if verbose:
                    print(f"  -> {classification['type']} (confidence: {classification['confidence']:.3f})")
                    print(f"     Reads: {classification['host_reads']} host, "
                          f"{classification['viral_reads']} viral, "
                          f"{classification['chimeric_reads']} chimeric")
        
        # Close BAM files
        host_bam.close()
        combined_bam.close()
        
        # Write outputs
        self._write_classified_bed(classifications, f"{output_prefix}_classified.bed")
        self._write_summary_report(classifications, f"{output_prefix}_summary.txt")
        
        if verbose:
            self._print_summary_statistics(classifications)
        
        return classifications
    
    def _classify_region(self, chrom: str, start: int, end: int, name: str,
                        host_bam: pysam.AlignmentFile, 
                        combined_bam: pysam.AlignmentFile,
                        min_mapq: int) -> Dict:
        """Classify a single region"""
        
        # Analyze read composition
        read_stats = self._analyze_reads(
            chrom, start, end, host_bam, combined_bam, min_mapq
        )
        
        # Detect integration junctions
        junction_info = self._detect_junctions(
            chrom, start, end, combined_bam, min_mapq
        )
        
        # Classify based on evidence
        region_type, confidence = self._determine_type(read_stats, junction_info)
        
        # Calculate additional metrics
        gc_content = self._calculate_gc(chrom, start, end)
        
        return {
            'name': name,
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': end - start,
            'type': region_type,
            'confidence': confidence,
            'total_reads': read_stats['total'],
            'host_reads': read_stats['host_only'],
            'viral_reads': read_stats['viral_only'],
            'chimeric_reads': read_stats['chimeric'],
            'unmapped_reads': read_stats['unmapped'],
            'has_junction': junction_info['has_junction'],
            'junction_reads': junction_info['junction_reads'],
            'gc_content': gc_content,
            'viral_species': read_stats.get('viral_species', [])
        }
    
    def _analyze_reads(self, chrom: str, start: int, end: int,
                      host_bam: pysam.AlignmentFile,
                      combined_bam: pysam.AlignmentFile,
                      min_mapq: int) -> Dict:
        """Analyze read composition in region"""
        
        stats = {
            'total': 0,
            'host_only': 0,
            'viral_only': 0,
            'chimeric': 0,
            'unmapped': 0,
            'viral_species': set()
        }
        
        # Get unique reads from region
        read_ids = set()
        for read in combined_bam.fetch(chrom, start, end):
            if read.mapping_quality >= min_mapq and not read.is_secondary:
                read_ids.add(read.query_name)
        
        # Classify each read
        for read_id in read_ids:
            stats['total'] += 1
            
            # Check if maps to host-only reference
            maps_to_host = self._read_in_host_bam(read_id, host_bam, min_mapq)
            
            # Check alignments in combined reference
            alignments = self._get_read_alignments(read_id, combined_bam, min_mapq)
            
            has_host = any(not ref.startswith(self.viral_prefix) 
                          for ref in alignments['refs'])
            has_viral = any(ref.startswith(self.viral_prefix) 
                           for ref in alignments['refs'])
            
            # Classify read
            if has_host and has_viral:
                stats['chimeric'] += 1
            elif has_viral and not maps_to_host:
                stats['viral_only'] += 1
                # Track viral species
                for ref in alignments['refs']:
                    if ref.startswith(self.viral_prefix):
                        species = ref.replace(self.viral_prefix, '').strip('_')
                        stats['viral_species'].add(species)
            elif has_host or maps_to_host:
                stats['host_only'] += 1
            else:
                stats['unmapped'] += 1
        
        stats['viral_species'] = list(stats['viral_species'])
        return stats
    
    def _detect_junctions(self, chrom: str, start: int, end: int,
                         combined_bam: pysam.AlignmentFile,
                         min_mapq: int) -> Dict:
        """Detect host-viral integration junctions"""
        
        junction_reads = 0
        breakpoints = []
        
        for read in combined_bam.fetch(chrom, start, end):
            if read.mapping_quality < min_mapq or read.is_secondary:
                continue
            
            # Check for supplementary alignments to viral
            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                for sa in sa_tag.rstrip(';').split(';'):
                    sa_ref = sa.split(',')[0]
                    if sa_ref.startswith(self.viral_prefix):
                        junction_reads += 1
                        breakpoints.append(read.reference_start)
            
            # Check for significant soft clipping
            if read.cigartuples:
                if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 20:
                    junction_reads += 1
                if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 20:
                    junction_reads += 1
        
        return {
            'has_junction': junction_reads > 0,
            'junction_reads': junction_reads,
            'breakpoints': breakpoints
        }
    
    def _determine_type(self, read_stats: Dict, junction_info: Dict) -> Tuple[str, float]:
        """Determine region type based on evidence"""
        
        total = read_stats['total']
        if total == 0:
            return 'unknown', 0.0
        
        # Calculate proportions
        host_prop = read_stats['host_only'] / total
        viral_prop = read_stats['viral_only'] / total
        chimeric_prop = read_stats['chimeric'] / total
        
        # Classification logic
        if junction_info['has_junction'] and chimeric_prop > 0.1:
            # Integration site
            confidence = min(0.95, chimeric_prop + 0.3)
            return 'integration_site', confidence
        
        elif viral_prop > 0.8:
            # Pure viral
            return 'viral_only', viral_prop
        
        elif host_prop > 0.8:
            # Pure host
            return 'host_only', host_prop
        
        elif chimeric_prop > 0.3 or (viral_prop > 0.2 and host_prop > 0.2):
            # Chimeric
            confidence = max(chimeric_prop, (viral_prop + host_prop) / 2)
            return 'chimeric', confidence
        
        else:
            # Mixed/unclear
            return 'mixed', 0.5
    
    def _read_in_host_bam(self, read_id: str, host_bam: pysam.AlignmentFile,
                         min_mapq: int) -> bool:
        """Check if read exists in host-only BAM"""
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
        alignments = {'refs': [], 'positions': []}
        
        for read in bam.fetch(until_eof=True):
            if read.query_name == read_id and read.mapping_quality >= min_mapq:
                alignments['refs'].append(read.reference_name)
                alignments['positions'].append(read.reference_start)
        
        return alignments
    
    def _calculate_gc(self, chrom: str, start: int, end: int) -> float:
        """Calculate GC content"""
        try:
            if chrom in self.combined_ref.references:
                seq = self.combined_ref.fetch(chrom, start, end).upper()
            else:
                seq = self.host_ref.fetch(chrom, start, end).upper()
            
            gc_count = seq.count('G') + seq.count('C')
            return gc_count / len(seq) if len(seq) > 0 else 0.0
        except:
            return 0.0
    
    def _write_classified_bed(self, classifications: List[Dict], output_file: str):
        """Write classified regions to BED file"""
        with open(output_file, 'w') as f:
            # Header
            f.write("#chrom\tstart\tend\tname\tscore\tstrand\tclassification\t")
            f.write("confidence\thost_reads\tviral_reads\tchimeric_reads\t")
            f.write("has_junction\tviral_species\n")
            
            for c in classifications:
                # Convert confidence to BED score (0-1000)
                score = int(c['confidence'] * 1000)
                
                # Format viral species
                viral_str = ','.join(c['viral_species']) if c['viral_species'] else 'none'
                
                # Write line
                f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t")
                f.write(f"{c['name']}\t{score}\t.\t{c['type']}\t")
                f.write(f"{c['confidence']:.3f}\t{c['host_reads']}\t")
                f.write(f"{c['viral_reads']}\t{c['chimeric_reads']}\t")
                f.write(f"{c['has_junction']}\t{viral_str}\n")
    
    def _write_summary_report(self, classifications: List[Dict], output_file: str):
        """Write summary report"""
        with open(output_file, 'w') as f:
            f.write("CircONTrack Classification Summary\n")
            f.write("=" * 50 + "\n\n")
            
            # Overall statistics
            total = len(classifications)
            f.write(f"Total regions classified: {total}\n\n")
            
            # Count by type
            type_counts = defaultdict(int)
            for c in classifications:
                type_counts[c['type']] += 1
            
            f.write("Classification breakdown:\n")
            for region_type in ['host_only', 'viral_only', 'chimeric', 'integration_site', 'mixed', 'unknown']:
                count = type_counts.get(region_type, 0)
                percent = (count / total * 100) if total > 0 else 0
                f.write(f"  {region_type:20s}: {count:4d} ({percent:5.1f}%)\n")
            
            # Integration sites
            integration_sites = [c for c in classifications if c['type'] == 'integration_site']
            if integration_sites:
                f.write(f"\n{len(integration_sites)} Integration sites detected:\n")
                for site in integration_sites:
                    f.write(f"  - {site['name']}: {site['chrom']}:{site['start']}-{site['end']}")
                    f.write(f" ({site['junction_reads']} junction reads)\n")
            
            # Viral species detected
            all_species = set()
            for c in classifications:
                all_species.update(c.get('viral_species', []))
            
            if all_species:
                f.write(f"\nViral species detected: {', '.join(sorted(all_species))}\n")
            
            # High-confidence regions
            high_conf = [c for c in classifications if c['confidence'] >= 0.8]
            f.write(f"\nHigh-confidence regions (≥0.8): {len(high_conf)}\n")
    
    def _print_summary_statistics(self, classifications: List[Dict]):
        """Print summary to console"""
        print("\n" + "=" * 50)
        print("Classification Summary")
        print("=" * 50)
        
        total = len(classifications)
        print(f"Total regions: {total}")
        
        # Count by type
        type_counts = defaultdict(int)
        for c in classifications:
            type_counts[c['type']] += 1
        
        print("\nBreakdown by type:")
        for region_type, count in sorted(type_counts.items()):
            percent = (count / total * 100) if total > 0 else 0
            print(f"  {region_type:20s}: {count:3d} ({percent:5.1f}%)")
        
        # Integration sites
        integration_count = type_counts.get('integration_site', 0)
        if integration_count > 0:
            print(f"\n⚠️  {integration_count} potential viral integration site(s) detected!")
        
        print("=" * 50)


def main():
    """Main entry point for classify subcommand"""
    parser = argparse.ArgumentParser(
        description='Classify CircONTrack-detected circular DNA regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic classification
  circontrack classify circdna.bed host.fa combined.fa host.bam combined.bam -o results
  
  # With custom viral prefix
  circontrack classify circdna.bed host.fa combined.fa host.bam combined.bam \\
    --viral-prefix "virus_" -o results
  
  # Quiet mode with higher quality threshold
  circontrack classify circdna.bed host.fa combined.fa host.bam combined.bam \\
    --min-mapq 30 --quiet -o results
        """
    )
    
    # Required arguments
    parser.add_argument('bed_file', 
                       help='BED file from CircONTrack detection')
    parser.add_argument('host_ref', 
                       help='Host-only reference FASTA (e.g., mm10)')
    parser.add_argument('combined_ref', 
                       help='Combined host+viral reference FASTA')
    parser.add_argument('host_bam', 
                       help='BAM file aligned to host-only reference')
    parser.add_argument('combined_bam', 
                       help='BAM file aligned to combined reference')
    
    # Optional arguments
    parser.add_argument('-o', '--output', default='classified',
                       help='Output prefix (default: classified)')
    parser.add_argument('--viral-prefix', default='viral',
                       help='Prefix for viral contigs in combined reference (default: viral)')
    parser.add_argument('--min-mapq', type=int, default=20,
                       help='Minimum mapping quality (default: 20)')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress messages')
    
    args = parser.parse_args()
    
    # Validate input files exist
    for file_path, name in [(args.bed_file, 'BED file'),
                            (args.host_ref, 'Host reference'),
                            (args.combined_ref, 'Combined reference'),
                            (args.host_bam, 'Host BAM'),
                            (args.combined_bam, 'Combined BAM')]:
        if not os.path.exists(file_path):
            print(f"Error: {name} not found: {file_path}")
            sys.exit(1)
    
    # Initialize classifier
    print("Initializing classifier...")
    classifier = CircularDNAClassifier(
        args.host_ref,
        args.combined_ref,
        args.viral_prefix
    )
    
    # Run classification
    print(f"Classifying regions from {args.bed_file}...")
    classifications = classifier.classify_bed_file(
        args.bed_file,
        args.host_bam,
        args.combined_bam,
        args.output,
        args.min_mapq,
        verbose=not args.quiet
    )
    
    print(f"\nResults written to:")
    print(f"  - {args.output}_classified.bed")
    print(f"  - {args.output}_summary.txt")


if __name__ == "__main__":
    main()