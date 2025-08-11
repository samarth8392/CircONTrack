#!/usr/bin/env python3
"""
CircONTrack Classify Module - Updated for RefSeq viral sequences
Handles real viral contig names (NC_*, NR_*) and extracts species from FASTA headers
"""

import pysam
import numpy as np
import re
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Set
import argparse
import sys
import os

class CircularDNAClassifier:
    """
    Classify circular DNA using combined reference BAM
    Updated to handle RefSeq viral sequences (NC_*, NR_*, etc.)
    """
    
    def __init__(self, combined_ref_path: str, 
                 viral_patterns: List[str] = None,
                 parse_viral_names: bool = True):
        """
        Initialize with combined reference
        
        Args:
            combined_ref_path: Combined host+viral reference FASTA
            viral_patterns: List of patterns to identify viral contigs 
                           (default: ['NC_', 'NR_', 'NZ_', 'viral'])
            parse_viral_names: Extract virus names from FASTA headers
        """
        self.combined_ref = pysam.FastaFile(combined_ref_path)
        
        # Default patterns for viral contigs (RefSeq accessions)
        if viral_patterns is None:
            self.viral_patterns = ['NC_', 'NR_', 'NZ_', 'viral', 'virus']
        else:
            self.viral_patterns = viral_patterns
        
        # Build contig classification and name mapping
        self.viral_contigs = set()
        self.host_contigs = set()
        self.viral_contig_to_species = {}  # Map contig ID to species name
        
        self._classify_contigs(combined_ref_path, parse_viral_names)
        
        print(f"Reference contains {len(self.host_contigs)} host and {len(self.viral_contigs)} viral contigs")
        
        if self.viral_contigs and len(self.viral_contigs) <= 20:
            print("\nViral contigs detected:")
            for contig in sorted(self.viral_contigs)[:20]:
                species = self.viral_contig_to_species.get(contig, 'Unknown')
                print(f"  {contig}: {species}")
            if len(self.viral_contigs) > 20:
                print(f"  ... and {len(self.viral_contigs) - 20} more")
    
    def _classify_contigs(self, fasta_path: str, parse_names: bool):
        """
        Classify contigs as host or viral and extract species names
        
        Parses FASTA headers like:
        >NC_029549 High Plains wheat mosaic virus isolate Nebraska segment RNA 2, complete sequence
        >chr1
        """
        print("Parsing reference FASTA to identify viral contigs...")
        
        with open(fasta_path, 'r') as f:
            current_contig = None
            
            for line in f:
                if line.startswith('>'):
                    # Parse header
                    header = line[1:].strip()
                    
                    # Extract contig ID (first word)
                    contig_id = header.split()[0]
                    
                    # Check if viral based on patterns
                    is_viral = any(contig_id.startswith(pattern) for pattern in self.viral_patterns)
                    
                    if is_viral:
                        self.viral_contigs.add(contig_id)
                        
                        if parse_names and ' ' in header:
                            # Extract species name from header
                            # Remove contig ID to get description
                            description = header[len(contig_id):].strip()
                            
                            # Clean up the species name
                            species_name = self._extract_species_name(description)
                            self.viral_contig_to_species[contig_id] = species_name
                        else:
                            # Use contig ID as species name
                            self.viral_contig_to_species[contig_id] = contig_id
                    else:
                        # Host contig (chr*, scaffold*, etc.)
                        self.host_contigs.add(contig_id)
    
    def _extract_species_name(self, description: str) -> str:
        """
        Extract clean species name from FASTA description
        
        Examples:
        "High Plains wheat mosaic virus isolate Nebraska segment RNA 2, complete sequence"
        → "High Plains wheat mosaic virus"
        
        "Human papillomavirus type 16, complete genome"
        → "Human papillomavirus type 16"
        """
        # Remove common suffixes
        clean_name = description
        
        # Remove segment information
        clean_name = re.sub(r'(segment|segment\s+\w+|segment\s+RNA\s+\d+)', '', clean_name, flags=re.IGNORECASE)
        
        # Remove sequence descriptors
        clean_name = re.sub(r',?\s*(complete\s+sequence|complete\s+genome|partial\s+sequence|genome)', '', 
                           clean_name, flags=re.IGNORECASE)
        
        # Remove isolate/strain info (optional - keep if you want strain-level detail)
        clean_name = re.sub(r'(isolate|strain|clone)\s+[\w\-\.]+', '', clean_name, flags=re.IGNORECASE)
        
        # Clean up whitespace
        clean_name = ' '.join(clean_name.split())
        
        # If name is too long, truncate intelligently
        if len(clean_name) > 50:
            # Try to keep just virus name
            if 'virus' in clean_name.lower():
                parts = clean_name.split('virus')
                clean_name = parts[0] + 'virus'
            else:
                clean_name = clean_name[:50] + '...'
        
        return clean_name.strip() or description[:50]
    
    def _is_viral_contig(self, contig_name: str) -> bool:
        """Check if a contig is viral based on name patterns or known viral contigs"""
        # First check our pre-classified set
        if contig_name in self.viral_contigs:
            return True
        
        # Also check patterns in case of SA tags with contigs not in main reference
        return any(contig_name.startswith(pattern) for pattern in self.viral_patterns)
    
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
            'viral_contigs': set(),  # Track which viral contigs are present
            'viral_species': set(),   # Track virus species names
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
            read_class, viral_contigs = self._classify_read(read)
            
            if read_class == 'host_only':
                stats['host_only'] += 1
            elif read_class == 'viral_only':
                stats['viral_only'] += 1
                stats['viral_contigs'].update(viral_contigs)
            elif read_class == 'chimeric':
                stats['chimeric'] += 1
                stats['viral_contigs'].update(viral_contigs)
            elif read_class == 'integration':
                stats['integration_junction'] += 1
                stats['chimeric'] += 1  # Count as chimeric too
                stats['viral_contigs'].update(viral_contigs)
                
                # Record junction position
                if read.reference_start >= start and read.reference_start <= end:
                    stats['junction_positions'].append(read.reference_start)
        
        # Convert viral contigs to species names
        for contig in stats['viral_contigs']:
            species = self.viral_contig_to_species.get(contig, contig)
            stats['viral_species'].add(species)
        
        # Determine region type
        region_type, confidence = self._determine_region_type(stats)
        
        stats['type'] = region_type
        stats['confidence'] = confidence
        stats['viral_species'] = list(stats['viral_species'])
        stats['viral_contigs'] = list(stats['viral_contigs'])
        
        if verbose and stats['viral_species']:
            print(f"  Viral species detected: {', '.join(stats['viral_species'][:3])}")
            if len(stats['viral_species']) > 3:
                print(f"    ... and {len(stats['viral_species']) - 3} more")
        
        return stats
    
    def _classify_read(self, read) -> Tuple[str, Set[str]]:
        """
        Classify a single read based on its alignments
        
        Returns:
            (classification, set_of_viral_contigs)
        """
        viral_contigs = set()
        
        # Check primary alignment
        primary_is_viral = self._is_viral_contig(read.reference_name)
        if primary_is_viral:
            viral_contigs.add(read.reference_name)
        
        # Check for supplementary alignments
        has_host_supp = False
        has_viral_supp = False
        
        if read.has_tag('SA'):
            sa_tag = read.get_tag('SA')
            for sa_entry in sa_tag.rstrip(';').split(';'):
                if not sa_entry:
                    continue
                    
                sa_fields = sa_entry.split(',')
                if len(sa_fields) < 6:
                    continue
                    
                sa_ref = sa_fields[0]
                
                if self._is_viral_contig(sa_ref):
                    has_viral_supp = True
                    viral_contigs.add(sa_ref)
                else:
                    has_host_supp = True
        
        # Check for large soft clips (potential integration signature)
        has_large_clips = False
        if read.cigartuples:
            # Check first and last operations
            if len(read.cigartuples) > 0:
                if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 30:
                    has_large_clips = True
                if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 30:
                    has_large_clips = True
        
        # Classify based on alignment pattern
        if primary_is_viral:
            if has_host_supp:
                return ('integration', viral_contigs)
            else:
                return ('viral_only', viral_contigs)
        else:  # Primary is host
            if has_viral_supp:
                if has_large_clips:
                    return ('integration', viral_contigs)
                else:
                    return ('chimeric', viral_contigs)
            else:
                return ('host_only', viral_contigs)
    
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
            f.write("viral_species\tviral_contigs\n")
            
            for c in classifications:
                score = int(c['confidence'] * 1000)
                
                # Format viral species (clean names)
                viral_species_str = '|'.join(c['viral_species'][:3]) if c['viral_species'] else 'none'
                if len(c['viral_species']) > 3:
                    viral_species_str += f'_and_{len(c["viral_species"])-3}_more'
                
                # Format viral contigs (accessions)
                viral_contigs_str = ','.join(c['viral_contigs']) if c['viral_contigs'] else 'none'
                
                f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t")
                f.write(f"{c['name']}\t{score}\t.\t{c['type']}\t")
                f.write(f"{c['confidence']:.3f}\t{c['total_reads']}\t")
                f.write(f"{c['host_only']}\t{c['viral_only']}\t")
                f.write(f"{c['chimeric']}\t{c['integration_junction']}\t")
                f.write(f"{viral_species_str}\t{viral_contigs_str}\n")
        
        # Write detailed summary
        with open(f"{output_prefix}_summary.txt", 'w') as f:
            f.write("CircONTrack Classification Summary\n")
            f.write("=" * 80 + "\n\n")
            
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
            
            # Viral species summary
            all_species = defaultdict(int)
            all_contigs = defaultdict(int)
            
            for c in classifications:
                for species in c.get('viral_species', []):
                    all_species[species] += 1
                for contig in c.get('viral_contigs', []):
                    all_contigs[contig] += 1
            
            if all_species:
                f.write(f"\n{len(all_species)} viral species detected:\n")
                # Sort by frequency
                for species, count in sorted(all_species.items(), key=lambda x: -x[1])[:20]:
                    f.write(f"  {species}: {count} regions\n")
                if len(all_species) > 20:
                    f.write(f"  ... and {len(all_species) - 20} more species\n")
                
                f.write(f"\nViral contigs (RefSeq accessions):\n")
                for contig, count in sorted(all_contigs.items(), key=lambda x: -x[1])[:10]:
                    species = self.viral_contig_to_species.get(contig, 'Unknown')
                    f.write(f"  {contig}: {count} regions ({species})\n")
            
            # Integration sites
            integration_sites = [c for c in classifications if c['type'] == 'integration_site']
            if integration_sites:
                f.write(f"\n{len(integration_sites)} potential integration sites:\n")
                for site in integration_sites[:20]:
                    f.write(f"  {site['name']}: {site['chrom']}:{site['start']}-{site['end']}")
                    if site['viral_species']:
                        species_str = ', '.join(site['viral_species'][:2])
                        if len(site['viral_species']) > 2:
                            species_str += f" and {len(site['viral_species'])-2} more"
                        f.write(f" ({species_str})")
                    f.write(f" [{site['integration_junction']} junction reads]\n")
                if len(integration_sites) > 20:
                    f.write(f"  ... and {len(integration_sites) - 20} more sites\n")
    
    def _print_summary(self, classifications: List[Dict]):
        """Print summary to console"""
        print("\n" + "=" * 80)
        print("CLASSIFICATION SUMMARY")
        print("=" * 80)
        
        total = len(classifications)
        type_counts = defaultdict(int)
        for c in classifications:
            type_counts[c['type']] += 1
        
        for region_type, count in sorted(type_counts.items()):
            pct = count / total * 100
            print(f"{region_type:20s}: {count:3d} ({pct:5.1f}%)")
        
        # Collect all viral species
        all_species = set()
        for c in classifications:
            all_species.update(c.get('viral_species', []))
        
        if all_species:
            print(f"\n🦠 {len(all_species)} different viral species detected")
            # Show top 5
            species_counts = defaultdict(int)
            for c in classifications:
                for species in c.get('viral_species', []):
                    species_counts[species] += 1
            
            print("Top viral species:")
            for species, count in sorted(species_counts.items(), key=lambda x: -x[1])[:5]:
                print(f"  - {species}: {count} regions")
        
        # Highlight key findings
        if type_counts.get('integration_site', 0) > 0:
            print(f"\n⚠️  {type_counts['integration_site']} potential viral integration site(s) detected!")
        
        if type_counts.get('viral_episome', 0) > 0:
            print(f"🔬 {type_counts['viral_episome']} viral episome(s) detected")
        
        print("=" * 80)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Classify CircONTrack circular DNA regions (RefSeq-aware version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This version automatically detects RefSeq viral sequences (NC_*, NR_* accessions)
and extracts virus species names from FASTA headers.

Examples:
  # Basic usage with RefSeq viruses
  circontrack-classify circdna.bed combined_ref.fa combined.bam -o results
  
  # Custom patterns for viral contigs
  circontrack-classify circdna.bed combined_ref.fa combined.bam \\
    --viral-patterns NC_ NR_ viral_ phage_ -o results
  
  # Use contig names only (don't parse species from headers)
  circontrack-classify circdna.bed combined_ref.fa combined.bam \\
    --no-parse-names -o results
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
    parser.add_argument('--viral-patterns', nargs='+',
                       default=['NC_', 'NR_', 'NZ_', 'viral', 'virus'],
                       help='Patterns to identify viral contigs (default: NC_ NR_ NZ_ viral virus)')
    parser.add_argument('--no-parse-names', action='store_true',
                       help='Do not parse virus names from FASTA headers')
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
    classifier = CircularDNAClassifier(
        args.combined_ref,
        viral_patterns=args.viral_patterns,
        parse_viral_names=not args.no_parse_names
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