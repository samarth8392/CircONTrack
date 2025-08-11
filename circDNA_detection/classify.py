#!/usr/bin/env python3
"""
Simplified CircONTrack Classify Module
Only requires combined reference BAM since viral genomes are separate contigs
"""

import pysam
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
import argparse
import sys
import os

class SimplifiedCircularDNAClassifier:
    """
    Classify circular DNA using only combined reference BAM
    Works because viral genomes are separate contigs, not integrated into host chromosomes
    """
    
    def __init__(self, combined_ref_path: str, viral_contig_prefix: str = "viral"):
        """
        Initialize with just the combined reference
        
        Args:
            combined_ref_path: Combined host+viral reference FASTA
            viral_contig_prefix: Prefix identifying viral contigs
        """
        self.combined_ref = pysam.FastaFile(combined_ref_path)
        self.viral_prefix = viral_contig_prefix
        
        # Identify viral vs host contigs
        self.viral_contigs = set()
        self.host_contigs = set()
        
        for ref in self.combined_ref.references:
            if ref.startswith(viral_contig_prefix):
                self.viral_contigs.add(ref)
            else:
                self.host_contigs.add(ref)
        
        print(f"Reference contains {len(self.host_contigs)} host and {len(self.viral_contigs)} viral contigs")
        if self.viral_contigs:
            print(f"Viral contigs: {', '.join(sorted(self.viral_contigs))}")
    
    def classify_bed_file(self, bed_file: str, combined_bam_path: str, 
                          output_prefix: str, min_mapq: int = 20, 
                          verbose: bool = True):
        """
        Classify all regions from a BED file
        
        Args:
            bed_file: CircONTrack output BED file
            combined_bam_path: BAM aligned to combined host+viral reference
            output_prefix: Prefix for output files
            min_mapq: Minimum mapping quality
            verbose: Print progress
        """
        combined_bam = pysam.AlignmentFile(combined_bam_path, "rb")
        classifications = []
        
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
                    print(f"\nClassifying {name}: {chrom}:{start}-{end}")
                
                # Classify region
                result = self._classify_region(
                    chrom, start, end, name, combined_bam, min_mapq, verbose
                )
                classifications.append(result)
        
        combined_bam.close()
        
        # Write outputs
        self._write_results(classifications, output_prefix)
        
        if verbose:
            self._print_summary(classifications)
        
        return classifications
    
    def _classify_region(self, chrom: str, start: int, end: int, name: str,
                        bam: pysam.AlignmentFile, min_mapq: int, 
                        verbose: bool) -> Dict:
        """Classify a single region based on read alignments"""
        
        stats = {
            'name': name,
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': end - start,
            'total_reads': 0,
            'host_only': 0,
            'viral_only': 0,
            'chimeric': 0,
            'integration_junction': 0,
            'viral_species': set(),
            'junction_positions': []
        }
        
        # Analyze each read in the region
        read_ids_seen = set()
        
        for read in bam.fetch(chrom, start, end):
            if read.is_secondary or read.mapping_quality < min_mapq:
                continue
            
            # Avoid counting same read multiple times
            if read.query_name in read_ids_seen:
                continue
            read_ids_seen.add(read.query_name)
            
            stats['total_reads'] += 1
            
            # Classify this read
            read_class, viral_refs = self._classify_read(read)
            
            if read_class == 'host_only':
                stats['host_only'] += 1
            elif read_class == 'viral_only':
                stats['viral_only'] += 1
                stats['viral_species'].update(viral_refs)
            elif read_class == 'chimeric':
                stats['chimeric'] += 1
                stats['viral_species'].update(viral_refs)
            elif read_class == 'integration':
                stats['integration_junction'] += 1
                stats['chimeric'] += 1  # Count as chimeric too
                stats['viral_species'].update(viral_refs)
                
                # Record junction position
                if read.reference_start >= start and read.reference_start <= end:
                    stats['junction_positions'].append(read.reference_start)
        
        # Determine region type
        region_type, confidence = self._determine_region_type(stats)
        
        stats['type'] = region_type
        stats['confidence'] = confidence
        stats['viral_species'] = list(stats['viral_species'])
        
        if verbose and stats['viral_species']:
            print(f"  Viral species detected: {', '.join(stats['viral_species'])}")
        
        return stats
    
    def _classify_read(self, read) -> Tuple[str, set]:
        """
        Classify a single read based on its alignments
        
        Returns:
            (classification, set_of_viral_contigs)
        """
        viral_refs = set()
        
        # Check primary alignment
        primary_is_viral = read.reference_name in self.viral_contigs
        if primary_is_viral:
            viral_refs.add(read.reference_name)
        
        # Check for supplementary alignments
        has_host_supp = False
        has_viral_supp = False
        
        if read.has_tag('SA'):
            sa_tag = read.get_tag('SA')
            for sa_entry in sa_tag.rstrip(';').split(';'):
                if not sa_entry:
                    continue
                    
                sa_fields = sa_entry.split(',')
                sa_ref = sa_fields[0]
                
                if sa_ref in self.viral_contigs:
                    has_viral_supp = True
                    viral_refs.add(sa_ref)
                elif sa_ref in self.host_contigs:
                    has_host_supp = True
        
        # Check for large soft clips (potential integration signature)
        has_large_clips = False
        if read.cigartuples:
            # Check first and last operations
            if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 30:
                has_large_clips = True
            if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 30:
                has_large_clips = True
        
        # Classify based on alignment pattern
        if primary_is_viral:
            if has_host_supp:
                return ('integration', viral_refs)
            else:
                return ('viral_only', viral_refs)
        else:  # Primary is host
            if has_viral_supp:
                if has_large_clips:
                    return ('integration', viral_refs)
                else:
                    return ('chimeric', viral_refs)
            else:
                return ('host_only', viral_refs)
    
    def _determine_region_type(self, stats: Dict) -> Tuple[str, float]:
        """Determine region type from read statistics"""
        
        total = stats['total_reads']
        if total == 0:
            return ('unknown', 0.0)
        
        # Calculate proportions
        host_prop = stats['host_only'] / total
        viral_prop = stats['viral_only'] / total
        chimeric_prop = stats['chimeric'] / total
        integration_prop = stats['integration_junction'] / total
        
        # Decision logic
        if integration_prop > 0.1 or (chimeric_prop > 0.2 and len(stats['junction_positions']) > 0):
            # Integration site
            confidence = min(0.95, integration_prop + chimeric_prop)
            return ('integration_site', confidence)
        
        elif viral_prop > 0.7:
            # Viral episome
            if viral_prop > 0.9:
                return ('viral_episome', viral_prop)
            else:
                return ('viral_dominant', viral_prop)
        
        elif host_prop > 0.8:
            # Host eccDNA
            return ('host_eccDNA', host_prop)
        
        elif chimeric_prop > 0.3:
            # Chimeric/complex
            return ('chimeric', chimeric_prop)
        
        else:
            # Mixed/unclear
            confidence = max(host_prop, viral_prop, chimeric_prop)
            return ('mixed', confidence)
    
    def _write_results(self, classifications: List[Dict], output_prefix: str):
        """Write classification results"""
        
        # Write enhanced BED
        with open(f"{output_prefix}_classified.bed", 'w') as f:
            f.write("#chrom\tstart\tend\tname\tscore\tstrand\ttype\tconfidence\t")
            f.write("total_reads\thost_only\tviral_only\tchimeric\tintegration_junction\t")
            f.write("viral_species\n")
            
            for c in classifications:
                score = int(c['confidence'] * 1000)
                viral_str = ','.join(c['viral_species']) if c['viral_species'] else 'none'
                
                f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t")
                f.write(f"{c['name']}\t{score}\t.\t{c['type']}\t")
                f.write(f"{c['confidence']:.3f}\t{c['total_reads']}\t")
                f.write(f"{c['host_only']}\t{c['viral_only']}\t")
                f.write(f"{c['chimeric']}\t{c['integration_junction']}\t")
                f.write(f"{viral_str}\n")
        
        # Write summary
        with open(f"{output_prefix}_summary.txt", 'w') as f:
            f.write("CircONTrack Classification Summary\n")
            f.write("=" * 60 + "\n\n")
            
            total = len(classifications)
            f.write(f"Total regions analyzed: {total}\n\n")
            
            # Type breakdown
            type_counts = defaultdict(int)
            for c in classifications:
                type_counts[c['type']] += 1
            
            f.write("Classification breakdown:\n")
            for region_type in ['host_eccDNA', 'viral_episome', 'viral_dominant', 
                              'integration_site', 'chimeric', 'mixed', 'unknown']:
                count = type_counts.get(region_type, 0)
                if count > 0:
                    pct = count / total * 100
                    f.write(f"  {region_type:20s}: {count:4d} ({pct:5.1f}%)\n")
            
            # Viral species
            all_species = set()
            for c in classifications:
                all_species.update(c.get('viral_species', []))
            
            if all_species:
                f.write(f"\nViral species detected:\n")
                for species in sorted(all_species):
                    # Count regions with this species
                    count = sum(1 for c in classifications if species in c.get('viral_species', []))
                    f.write(f"  {species}: {count} regions\n")
            
            # Integration sites
            integration_sites = [c for c in classifications if c['type'] == 'integration_site']
            if integration_sites:
                f.write(f"\n{len(integration_sites)} potential integration sites:\n")
                for site in integration_sites:
                    f.write(f"  {site['name']}: {site['chrom']}:{site['start']}-{site['end']}")
                    if site['viral_species']:
                        f.write(f" ({', '.join(site['viral_species'])})")
                    f.write(f" [{site['integration_junction']} junction reads]\n")
    
    def _print_summary(self, classifications: List[Dict]):
        """Print summary to console"""
        print("\n" + "=" * 60)
        print("CLASSIFICATION SUMMARY")
        print("=" * 60)
        
        total = len(classifications)
        type_counts = defaultdict(int)
        for c in classifications:
            type_counts[c['type']] += 1
        
        for region_type, count in sorted(type_counts.items()):
            pct = count / total * 100
            print(f"{region_type:20s}: {count:3d} ({pct:5.1f}%)")
        
        # Highlight key findings
        if type_counts.get('integration_site', 0) > 0:
            print(f"\n⚠️  {type_counts['integration_site']} potential viral integration site(s) detected!")
        
        if type_counts.get('viral_episome', 0) > 0:
            print(f"🦠 {type_counts['viral_episome']} viral episome(s) detected")
        
        print("=" * 60)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Classify CircONTrack circular DNA regions (simplified version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This simplified version only requires:
  - Combined host+viral reference (with viral as separate contigs)
  - BAM file aligned to the combined reference
  
It does NOT need separate host-only alignments.

Examples:
  # Basic usage
  circontrack-classify circdna.bed combined_ref.fa combined.bam -o results
  
  # With custom viral prefix
  circontrack-classify circdna.bed combined_ref.fa combined.bam \\
    --viral-prefix "virus_" -o results
        """
    )
    
    parser.add_argument('bed_file', 
                       help='BED file from CircONTrack detection')
    parser.add_argument('combined_ref', 
                       help='Combined host+viral reference FASTA')
    parser.add_argument('combined_bam', 
                       help='BAM aligned to combined reference')
    
    parser.add_argument('-o', '--output', default='classified',
                       help='Output prefix (default: classified)')
    parser.add_argument('--viral-prefix', default='viral',
                       help='Prefix for viral contigs (default: "viral")')
    parser.add_argument('--min-mapq', type=int, default=20,
                       help='Minimum mapping quality (default: 20)')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress messages')
    
    args = parser.parse_args()
    
    # Validate files
    for filepath, name in [(args.bed_file, 'BED file'),
                           (args.combined_ref, 'Reference'),
                           (args.combined_bam, 'BAM file')]:
        if not os.path.exists(filepath):
            print(f"Error: {name} not found: {filepath}")
            sys.exit(1)
    
    # Run classification
    print("Initializing classifier...")
    classifier = SimplifiedCircularDNAClassifier(
        args.combined_ref,
        args.viral_prefix
    )
    
    print(f"Classifying regions from {args.bed_file}...")
    classifications = classifier.classify_bed_file(
        args.bed_file,
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