#!/usr/bin/env python3
"""
CircONTrack Viral Episome Detection Module
Detects viral episomes through uniform coverage analysis rather than peak detection
Part of the CircONTrack suite for circular DNA detection and analysis
"""

import pysam
import numpy as np
import pandas as pd
import re
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Set
import argparse
import sys
import os
from pathlib import Path
from circDNA_detection import __version__

class ViralEpisomeDetector:
    """
    Detect viral episomes by analyzing uniform coverage across viral contigs
    
    Unlike host eccDNA which appear as localized coverage peaks, viral episomes
    typically show uniform high coverage across the entire viral genome.
    """
    
    def __init__(self, bam_file: str, reference_file: str,
                 viral_patterns: List[str] = None,
                 min_coverage: float = 10.0,
                 min_uniformity: float = 0.7,
                 min_reads: int = 50,
                 min_mapq: int = 20,
                 parse_viral_names: bool = True,
                 verbose: bool = True):
        """
        Initialize viral episome detector
        
        Args:
            bam_file: Path to BAM file aligned to combined host+viral reference
            reference_file: Path to combined reference FASTA
            viral_patterns: Patterns to identify viral contigs (default: NC_, NR_, NZ_, viral, virus)
            min_coverage: Minimum mean coverage to call episome (default: 10)
            min_uniformity: Minimum uniformity score (default: 0.7)
            min_reads: Minimum total reads (default: 50)
            min_mapq: Minimum mapping quality (default: 20)
            parse_viral_names: Extract species names from FASTA headers (default: True)
            verbose: Print progress messages (default: True)
        """
        
        self.bam_file = bam_file
        self.reference_file = reference_file
        self.min_coverage = min_coverage
        self.min_uniformity = min_uniformity
        self.min_reads = min_reads
        self.min_mapq = min_mapq
        self.verbose = verbose
        
        # Default viral patterns
        if viral_patterns is None:
            self.viral_patterns = ['NC_', 'NR_', 'NZ_', 'viral', 'virus']
        else:
            self.viral_patterns = viral_patterns
        
        # Open files
        self.bam = pysam.AlignmentFile(bam_file, "rb")
        self.reference = pysam.FastaFile(reference_file)
        
        # Identify viral contigs
        self.viral_contigs = {}  # contig_id -> species_name
        self._identify_viral_contigs(parse_viral_names)
        
        if self.verbose:
            print(f"CircONTrack Viral Episome Detector v{__version__}")
            print(f"Identified {len(self.viral_contigs)} viral contigs in reference")
            if len(self.viral_contigs) <= 20:
                print("\nViral contigs:")
                for contig, species in sorted(self.viral_contigs.items()):
                    print(f"  {contig}: {species}")
    
    def _identify_viral_contigs(self, parse_names: bool):
        """
        Identify viral contigs from reference FASTA
        Reuses logic from circontrack-classify
        """
        
        if self.verbose:
            print("Parsing reference to identify viral contigs...")
        
        with open(self.reference_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    contig_id = header.split()[0]
                    
                    # Check if viral
                    is_viral = any(contig_id.startswith(p) for p in self.viral_patterns)
                    
                    if is_viral:
                        if parse_names and ' ' in header:
                            description = header[len(contig_id):].strip()
                            species_name = self._extract_species_name(description)
                            self.viral_contigs[contig_id] = species_name
                        else:
                            self.viral_contigs[contig_id] = contig_id
    
    def _extract_species_name(self, description: str) -> str:
        """Extract clean species name from FASTA description"""
        
        clean_name = description
        
        # Remove segment information
        clean_name = re.sub(r'(segment|segment\s+\w+|segment\s+RNA\s+\d+)', '', 
                           clean_name, flags=re.IGNORECASE)
        
        # Remove sequence descriptors
        clean_name = re.sub(r',?\s*(complete\s+sequence|complete\s+genome|partial\s+sequence|genome)', 
                           '', clean_name, flags=re.IGNORECASE)
        
        # Remove isolate/strain info
        clean_name = re.sub(r'(isolate|strain|clone)\s+[\w\-\.]+', '', 
                           clean_name, flags=re.IGNORECASE)
        
        # Clean whitespace
        clean_name = ' '.join(clean_name.split())
        
        # Truncate if too long
        if len(clean_name) > 50:
            if 'virus' in clean_name.lower():
                parts = clean_name.split('virus')
                clean_name = parts[0] + 'virus'
            else:
                clean_name = clean_name[:50] + '...'
        
        return clean_name.strip() or description[:50]
    
    def analyze_viral_contig(self, contig: str) -> Dict:
        """
        Analyze a single viral contig for episome evidence
        
        Args:
            contig: Viral contig name
            
        Returns:
            Dictionary with analysis results
        """
        
        if self.verbose:
            print(f"\nAnalyzing {contig} ({self.viral_contigs[contig]})...")
        
        # Get contig length
        contig_length = self.reference.get_reference_length(contig)
        
        result = {
            'contig': contig,
            'species': self.viral_contigs[contig],
            'length': contig_length,
            'mean_coverage': 0.0,
            'median_coverage': 0.0,
            'coverage_std': 0.0,
            'uniformity_score': 0.0,
            'cv': 0.0,
            'total_reads': 0,
            'high_quality_reads': 0,
            'mean_mapq': 0.0,
            'junction_reads': 0,
            'circular_junctions': 0,
            'is_circular': False,
            'is_episome': False,
            'episome_type': 'none',
            'confidence': 0.0
        }
        
        # Calculate coverage across contig
        coverage_data = self._calculate_coverage(contig)
        
        if len(coverage_data['coverage_array']) == 0:
            result['episome_type'] = 'no_coverage'
            return result
        
        result['mean_coverage'] = coverage_data['mean_coverage']
        result['median_coverage'] = coverage_data['median_coverage']
        result['coverage_std'] = coverage_data['std_coverage']
        result['cv'] = coverage_data['cv']
        result['uniformity_score'] = coverage_data['uniformity_score']
        
        # Analyze reads
        read_stats = self._analyze_reads(contig)
        result['total_reads'] = read_stats['total_reads']
        result['high_quality_reads'] = read_stats['high_quality_reads']
        result['mean_mapq'] = read_stats['mean_mapq']
        
        # Check for circular junctions
        junction_stats = self._detect_circular_junctions(contig, contig_length)
        result['junction_reads'] = junction_stats['junction_reads']
        result['circular_junctions'] = junction_stats['circular_junctions']
        result['is_circular'] = junction_stats['is_circular']
        
        # Classify episome
        self._classify_episome(result)
        
        if self.verbose:
            print(f"  Coverage: {result['mean_coverage']:.1f}x (uniformity: {result['uniformity_score']:.2f})")
            print(f"  Reads: {result['total_reads']} total, {result['high_quality_reads']} high-quality")
            print(f"  Type: {result['episome_type']}")
            if result['is_circular']:
                print(f"  ✓ Circular episome ({result['circular_junctions']} junction reads)")
        
        return result
    
    def _calculate_coverage(self, contig: str, bin_size: int = 100) -> Dict:
        """Calculate coverage statistics across viral contig"""
        
        contig_length = self.reference.get_reference_length(contig)
        
        # Calculate per-position coverage
        coverage_array = np.zeros(contig_length, dtype=np.int32)
        
        try:
            for pileup_column in self.bam.pileup(contig, 0, contig_length, 
                                                  truncate=True, 
                                                  max_depth=100000):
                coverage_array[pileup_column.pos] = pileup_column.n
        except Exception as e:
            print(f"Warning: Error calculating coverage for {contig}: {e}")
            return {
                'coverage_array': np.array([]),
                'mean_coverage': 0,
                'median_coverage': 0,
                'std_coverage': 0,
                'cv': 0,
                'uniformity_score': 0
            }
        
        # Calculate statistics
        mean_cov = np.mean(coverage_array)
        median_cov = np.median(coverage_array)
        std_cov = np.std(coverage_array)
        
        # Coefficient of variation
        cv = std_cov / mean_cov if mean_cov > 0 else 0
        
        # Uniformity score (inverse of CV, scaled 0-1)
        uniformity_score = 1 / (1 + cv) if cv > 0 else 1.0
        
        return {
            'coverage_array': coverage_array,
            'mean_coverage': mean_cov,
            'median_coverage': median_cov,
            'std_coverage': std_cov,
            'cv': cv,
            'uniformity_score': uniformity_score
        }
    
    def _analyze_reads(self, contig: str) -> Dict:
        """Analyze read mapping characteristics"""
        
        stats = {
            'total_reads': 0,
            'high_quality_reads': 0,
            'mapq_scores': [],
            'mean_mapq': 0.0
        }
        
        seen_reads = set()
        
        try:
            for read in self.bam.fetch(contig):
                # Count unique reads only
                if read.query_name in seen_reads:
                    continue
                seen_reads.add(read.query_name)
                
                if read.is_unmapped or read.is_secondary:
                    continue
                
                stats['total_reads'] += 1
                stats['mapq_scores'].append(read.mapping_quality)
                
                if read.mapping_quality >= self.min_mapq:
                    stats['high_quality_reads'] += 1
        
        except Exception as e:
            print(f"Warning: Error analyzing reads for {contig}: {e}")
        
        if stats['mapq_scores']:
            stats['mean_mapq'] = np.mean(stats['mapq_scores'])
        
        return stats
    
    def _detect_circular_junctions(self, contig: str, contig_length: int, 
                                   window: int = 500) -> Dict:
        """
        Detect circular junctions by finding reads that span contig boundaries
        
        Circular episomes have reads mapping from end → beginning of same contig
        """
        
        junction_stats = {
            'junction_reads': 0,
            'circular_junctions': 0,
            'is_circular': False,
            'junction_positions': []
        }
        
        try:
            # Check reads near the end of contig
            for read in self.bam.fetch(contig, max(0, contig_length - window), contig_length):
                
                if read.is_unmapped or read.mapping_quality < self.min_mapq:
                    continue
                
                # Check for supplementary alignment (SA tag)
                if read.has_tag('SA'):
                    sa_tag = read.get_tag('SA')
                    
                    for sa_entry in sa_tag.rstrip(';').split(';'):
                        if not sa_entry:
                            continue
                        
                        sa_fields = sa_entry.split(',')
                        if len(sa_fields) < 2:
                            continue
                        
                        sa_contig = sa_fields[0]
                        sa_pos = int(sa_fields[1])
                        
                        # Check if SA maps to beginning of same contig
                        if sa_contig == contig and sa_pos < window:
                            junction_stats['circular_junctions'] += 1
                            junction_stats['junction_positions'].append((read.reference_start, sa_pos))
                
                # Check for large soft clips that might indicate junction
                if read.cigartuples:
                    # Check for soft clip at end of read
                    if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 30:
                        junction_stats['junction_reads'] += 1
            
            # Also check reads at beginning for complementary pattern
            for read in self.bam.fetch(contig, 0, window):
                
                if read.is_unmapped or read.mapping_quality < self.min_mapq:
                    continue
                
                if read.has_tag('SA'):
                    sa_tag = read.get_tag('SA')
                    
                    for sa_entry in sa_tag.rstrip(';').split(';'):
                        if not sa_entry:
                            continue
                        
                        sa_fields = sa_entry.split(',')
                        if len(sa_fields) < 2:
                            continue
                        
                        sa_contig = sa_fields[0]
                        sa_pos = int(sa_fields[1])
                        
                        # Check if SA maps to end of same contig
                        if sa_contig == contig and sa_pos > contig_length - window:
                            junction_stats['circular_junctions'] += 1
                            junction_stats['junction_positions'].append((read.reference_start, sa_pos))
                
                # Check for soft clip at beginning
                if read.cigartuples:
                    if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 30:
                        junction_stats['junction_reads'] += 1
        
        except Exception as e:
            print(f"Warning: Error detecting junctions for {contig}: {e}")
        
        # Determine if circular based on junction evidence
        junction_stats['is_circular'] = junction_stats['circular_junctions'] >= 2
        
        return junction_stats
    
    def _classify_episome(self, result: Dict):
        """
        Classify viral contig as episome based on all evidence
        
        Types:
        - circular_episome: High coverage, uniform, circular junctions
        - linear_episome: High coverage, uniform, no circular junctions
        - low_coverage: Below threshold
        - non_uniform: Coverage present but not uniform (integration/fragmented)
        - no_coverage: No reads mapping
        """
        
        # Check basic thresholds
        meets_coverage = result['mean_coverage'] >= self.min_coverage
        meets_uniformity = result['uniformity_score'] >= self.min_uniformity
        meets_reads = result['total_reads'] >= self.min_reads
        
        # Calculate confidence score
        confidence = 0.0
        
        if meets_coverage:
            confidence += 30
        if meets_uniformity:
            confidence += 30
        if meets_reads:
            confidence += 20
        if result['is_circular']:
            confidence += 20
        
        result['confidence'] = min(100, confidence)
        
        # Classification logic
        if meets_coverage and meets_uniformity and meets_reads:
            result['is_episome'] = True
            
            if result['is_circular']:
                result['episome_type'] = 'circular_episome'
            else:
                result['episome_type'] = 'linear_episome'
        
        elif result['mean_coverage'] < self.min_coverage:
            result['episome_type'] = 'low_coverage'
            result['is_episome'] = False
        
        elif result['uniformity_score'] < self.min_uniformity:
            result['episome_type'] = 'non_uniform'
            result['is_episome'] = False
        
        elif result['total_reads'] < self.min_reads:
            result['episome_type'] = 'insufficient_reads'
            result['is_episome'] = False
        
        else:
            result['episome_type'] = 'uncertain'
            result['is_episome'] = False
    
    def detect_all_episomes(self) -> pd.DataFrame:
        """
        Analyze all viral contigs and detect episomes
        
        Returns:
            DataFrame with results for all viral contigs
        """
        
        if self.verbose:
            print(f"\nAnalyzing {len(self.viral_contigs)} viral contigs...")
        
        results = []
        
        for contig in sorted(self.viral_contigs.keys()):
            result = self.analyze_viral_contig(contig)
            results.append(result)
        
        df = pd.DataFrame(results)
        return df
    
    def save_results(self, results_df: pd.DataFrame, output_prefix: str):
        """Save results in multiple formats"""
        
        # Full results TSV
        full_output = f"{output_prefix}_all_viral_contigs.tsv"
        results_df.to_csv(full_output, sep='\t', index=False)
        
        if self.verbose:
            print(f"\nSaved all viral contig analysis to: {full_output}")
        
        # BED file for episomes only (compatible with CircONTrack pipeline)
        episomes = results_df[results_df['is_episome']]
        
        if len(episomes) > 0:
            bed_output = f"{output_prefix}_viral_episomes.bed"
            
            with open(bed_output, 'w') as f:
                f.write("# CircONTrack Viral Episome Detection Results\n")
                f.write(f"# Detected {len(episomes)} viral episomes\n")
                f.write("#chr\tstart\tend\tname\tscore\tstrand\ttype\tcoverage\tuniformity\tcircular\n")
                
                for _, row in episomes.iterrows():
                    # BED format: full viral contig as single region
                    score = int(row['confidence'])
                    f.write(f"{row['contig']}\t0\t{row['length']}\t")
                    f.write(f"{row['species']}\t{score}\t.\t")
                    f.write(f"{row['episome_type']}\t{row['mean_coverage']:.1f}\t")
                    f.write(f"{row['uniformity_score']:.2f}\t")
                    f.write(f"{'yes' if row['is_circular'] else 'no'}\n")
            
            if self.verbose:
                print(f"Saved {len(episomes)} viral episomes to: {bed_output}")
        
        # Summary report
        summary_output = f"{output_prefix}_summary.txt"
        
        with open(summary_output, 'w') as f:
            f.write("CircONTrack Viral Episome Detection Summary\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Total viral contigs analyzed: {len(results_df)}\n")
            f.write(f"Viral episomes detected: {len(episomes)}\n\n")
            
            # Type breakdown
            type_counts = results_df['episome_type'].value_counts()
            f.write("Classification breakdown:\n")
            for episome_type, count in type_counts.items():
                f.write(f"  {episome_type}: {count}\n")
            
            # Episome details
            if len(episomes) > 0:
                f.write("\n" + "=" * 70 + "\n")
                f.write("DETECTED VIRAL EPISOMES\n")
                f.write("=" * 70 + "\n\n")
                
                for _, row in episomes.iterrows():
                    f.write(f"Virus: {row['species']}\n")
                    f.write(f"  Contig: {row['contig']}\n")
                    f.write(f"  Length: {row['length']:,} bp\n")
                    f.write(f"  Coverage: {row['mean_coverage']:.1f}x (uniformity: {row['uniformity_score']:.2f})\n")
                    f.write(f"  Reads: {row['total_reads']} total, {row['high_quality_reads']} high-quality\n")
                    f.write(f"  Type: {row['episome_type']}\n")
                    f.write(f"  Confidence: {row['confidence']:.0f}/100\n")
                    
                    if row['is_circular']:
                        f.write(f"  ✓ Circular episome confirmed ({row['circular_junctions']} junction reads)\n")
                    
                    f.write("\n")
            
            # Coverage statistics
            f.write("=" * 70 + "\n")
            f.write("COVERAGE STATISTICS\n")
            f.write("=" * 70 + "\n\n")
            
            covered = results_df[results_df['mean_coverage'] >= 1.0]
            f.write(f"Viral contigs with coverage: {len(covered)}/{len(results_df)}\n")
            
            if len(covered) > 0:
                f.write(f"  Mean coverage range: {covered['mean_coverage'].min():.1f}x - {covered['mean_coverage'].max():.1f}x\n")
                f.write(f"  Median coverage: {covered['mean_coverage'].median():.1f}x\n")
                
                f.write("\nTop 10 by coverage:\n")
                top_coverage = results_df.nlargest(10, 'mean_coverage')
                for _, row in top_coverage.iterrows():
                    if row['mean_coverage'] >= 1.0:
                        f.write(f"  {row['species'][:50]}: {row['mean_coverage']:.1f}x ")
                        f.write(f"({row['episome_type']})\n")
        
        if self.verbose:
            print(f"Saved summary to: {summary_output}")
        
        # Print console summary
        self._print_summary(results_df, episomes)
    
    def _print_summary(self, results_df: pd.DataFrame, episomes: pd.DataFrame):
        """Print summary to console"""
        
        print("\n" + "=" * 70)
        print("VIRAL EPISOME DETECTION SUMMARY")
        print("=" * 70)
        
        print(f"\nTotal viral contigs analyzed: {len(results_df)}")
        print(f"Viral episomes detected: {len(episomes)}")
        
        if len(episomes) > 0:
            circular = episomes[episomes['is_circular']]
            linear = episomes[~episomes['is_circular']]
            
            print(f"  - Circular episomes: {len(circular)}")
            print(f"  - Linear episomes: {len(linear)}")
            
            print("\nDetected episomes:")
            for _, row in episomes.iterrows():
                circ_marker = "🔵" if row['is_circular'] else "⚪"
                print(f"  {circ_marker} {row['species']}")
                print(f"     Coverage: {row['mean_coverage']:.1f}x, Confidence: {row['confidence']:.0f}/100")
        
        else:
            print("\nNo viral episomes detected with current thresholds.")
            
            # Suggest threshold adjustments
            low_cov = results_df[
                (results_df['mean_coverage'] < self.min_coverage) & 
                (results_df['mean_coverage'] >= 1.0)
            ]
            
            if len(low_cov) > 0:
                print(f"\nNote: {len(low_cov)} viral contigs have coverage below threshold ({self.min_coverage}x)")
                print("Consider adjusting --min-coverage if investigating low-abundance viruses")
        
        print("=" * 70)
    
    def close(self):
        """Close open file handles"""
        self.bam.close()
        self.reference.close()


def main():
    """Main entry point for circontrack-viral-episomes"""
    
    parser = argparse.ArgumentParser(
        description='CircONTrack Viral Episome Detection - Detect viral episomes via uniform coverage',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This module detects viral episomes by analyzing uniform coverage patterns across
viral contigs, rather than using peak detection designed for host genomes.

Examples:
  # Basic usage
  circontrack-viral-episomes combined.bam combined_ref.fa -o viral_results
  
  # Custom thresholds for low-abundance viruses
  circontrack-viral-episomes combined.bam combined_ref.fa \\
    --min-coverage 5 --min-uniformity 0.6 -o viral_results
  
  # Custom viral patterns
  circontrack-viral-episomes combined.bam combined_ref.fa \\
    --viral-patterns NC_ viral_ phage_ -o viral_results

Output files:
  - {prefix}_viral_episomes.bed: BED format (mergeable with host eccDNA peaks)
  - {prefix}_all_viral_contigs.tsv: Complete analysis of all viral contigs
  - {prefix}_summary.txt: Human-readable summary report
        """
    )
    
    parser.add_argument('bam_file', 
                       help='BAM file aligned to combined host+viral reference')
    parser.add_argument('reference_file', 
                       help='Combined host+viral reference FASTA')
    
    parser.add_argument('-o', '--output-prefix', default='viral_episomes',
                       help='Output file prefix (default: viral_episomes)')
    
    parser.add_argument('--viral-patterns', nargs='+',
                       default=['NC_', 'NR_', 'NZ_', 'viral', 'virus'],
                       help='Patterns to identify viral contigs (default: NC_ NR_ NZ_ viral virus)')
    
    parser.add_argument('--min-coverage', type=float, default=10.0,
                       help='Minimum mean coverage for episome (default: 10)')
    parser.add_argument('--min-uniformity', type=float, default=0.7,
                       help='Minimum uniformity score (default: 0.7)')
    parser.add_argument('--min-reads', type=int, default=50,
                       help='Minimum total reads (default: 50)')
    parser.add_argument('--min-mapq', type=int, default=20,
                       help='Minimum mapping quality (default: 20)')
    
    parser.add_argument('--no-parse-names', action='store_true',
                       help='Do not parse virus names from FASTA headers')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Suppress progress messages')
    
    parser.add_argument('--version', action='version', 
                       version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    # Validate input files
    for filepath, name in [(args.bam_file, 'BAM file'),
                           (args.reference_file, 'Reference FASTA')]:
        if not os.path.exists(filepath):
            print(f"Error: {name} not found: {filepath}")
            sys.exit(1)
    
    # Check BAM index
    bam_index = f"{args.bam_file}.bai"
    if not os.path.exists(bam_index):
        print(f"Error: BAM index not found: {bam_index}")
        print("Please index your BAM file with: samtools index <bam_file>")
        sys.exit(1)
    
    try:
        # Initialize detector
        detector = ViralEpisomeDetector(
            bam_file=args.bam_file,
            reference_file=args.reference_file,
            viral_patterns=args.viral_patterns,
            min_coverage=args.min_coverage,
            min_uniformity=args.min_uniformity,
            min_reads=args.min_reads,
            min_mapq=args.min_mapq,
            parse_viral_names=not args.no_parse_names,
            verbose=not args.quiet
        )
        
        # Detect episomes
        results = detector.detect_all_episomes()
        
        # Save results
        detector.save_results(results, args.output_prefix)
        
        # Cleanup
        detector.close()
        
        if not args.quiet:
            print(f"\nViral episome detection complete!")
            print(f"Output files: {args.output_prefix}_*")
    
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()