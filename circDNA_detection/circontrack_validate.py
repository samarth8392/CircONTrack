#!/usr/bin/env python3
"""
CircONTrack Peak Validation Module
Combines coverage peak calling with junction validation to identify true eccDNA
while filtering out repetitive element artifacts

This module is part of the CircONTrack suite for circular DNA detection and analysis.
"""

import pysam
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import logging
from pathlib import Path
import argparse
import sys
import os

# Add CircONTrack modules to path if running as part of the suite
try:
    from circDNA_detection.utils import setup_logging, validate_file
    from circDNA_detection.config import DEFAULT_PARAMS
except ImportError:
    # Fallback functions if not part of CircONTrack
    def setup_logging(verbose=False):
        import logging
        level = logging.DEBUG if verbose else logging.INFO
        logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)
    
    def validate_file(filepath):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        return True

__version__ = "1.0.0"
__author__ = "CircONTrack Development Team"


@dataclass
class CircularCandidate:
    """Enhanced circular DNA candidate with all evidence"""
    chromosome: str
    start: int
    end: int
    length: int
    
    # Coverage evidence
    fold_change: float = 0.0
    coverage: float = 0.0
    coverage_uniformity: float = 0.0
    pvalue: float = 1.0
    adjusted_pvalue: float = 1.0
    
    # Junction evidence  
    junction_reads: int = 0
    split_reads: int = 0
    discordant_pairs: int = 0
    boundary_clips: int = 0
    
    # Quality metrics
    mean_mapq: float = 0.0
    soft_clip_rate: float = 0.0
    position_diversity: float = 0.0
    
    # Classification
    confidence_score: float = 0.0
    classification: str = "unknown"
    evidence_summary: str = ""
    is_artifact: bool = False
    artifact_reasons: List[str] = None
    
    def __post_init__(self):
        if self.artifact_reasons is None:
            self.artifact_reasons = []


class CircONTrackValidator:
    """
    Integrated CircONTrack validation pipeline that combines:
    1. Coverage peak analysis (compatible with coverage_peaks.py output)
    2. Junction detection and validation
    3. Artifact filtering and classification
    4. Confidence scoring for eccDNA candidates
    
    This module validates peaks detected by CircONTrack's coverage_peaks.py module.
    """
    
    def __init__(self, bam_file, peak_file=None, reference_file=None, 
                 min_junction_support=2, max_soft_clip_rate=0.3,
                 min_mapq=20, max_fold_change=100, logger=None):
        """
        Initialize the CircONTrack validation pipeline
        
        Parameters:
        -----------
        bam_file : str
            Path to BAM file used for peak calling
        peak_file : str, optional
            Path to CircONTrack peak output file
        reference_file : str, optional
            Path to reference FASTA file
        min_junction_support : int
            Minimum junction reads required for eccDNA classification
        max_soft_clip_rate : float
            Maximum soft-clipping rate before flagging as artifact
        min_mapq : int
            Minimum mean mapping quality
        max_fold_change : float
            Maximum fold change before flagging as potential artifact
        logger : logging.Logger, optional
            Logger instance for output messages
        """
        
        self.bam_file = bam_file
        self.peak_file = peak_file
        self.reference_file = reference_file
        self.logger = logger or setup_logging()
        
        # Validate input files
        validate_file(self.bam_file)
        if self.peak_file:
            validate_file(self.peak_file)
        if self.reference_file:
            validate_file(self.reference_file)
        
        # Open BAM file
        try:
            self.bam = pysam.AlignmentFile(bam_file, 'rb')
            self.logger.info(f"Opened BAM file: {bam_file}")
        except Exception as e:
            self.logger.error(f"Error opening BAM file {bam_file}: {str(e)}")
            raise
        
        # Validation thresholds
        self.min_junction_support = min_junction_support
        self.max_soft_clip_rate = max_soft_clip_rate
        self.min_mapq = min_mapq
        self.max_fold_change = max_fold_change
        
        # Load peaks if provided
        self.peaks = None
        if peak_file:
            self.peaks = self.load_peaks(peak_file)
        
        self.logger.info("CircONTrack Validator initialized successfully")
    
    def load_peaks(self, peak_file):
        """
        Load peaks from CircONTrack coverage peak caller output
        
        Expected format matches CircONTrack's coverage_peaks.py output:
        chr, start, end, name, score, strand, coverage, fold_change, pvalue, adjusted_pvalue, read_count
        """
        
        self.logger.info(f"Loading peaks from {peak_file}")
        
        try:
            # Read the file, skipping comment lines that start with #
            peaks = pd.read_csv(
                peak_file, 
                sep='\t', 
                comment='#',
                low_memory=False
            )
            
            # Expected CircONTrack peak output columns
            expected_columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 
                              'coverage', 'fold_change', 'pvalue', 'adjusted_pvalue', 'read_count']
            
            # Handle header detection and column assignment
            if len(peaks.columns) != len(expected_columns):
                self.logger.warning(f"Expected {len(expected_columns)} columns, got {len(peaks.columns)}")
                
                # Try reading without skipping comments to see the actual header
                peaks_with_header = pd.read_csv(peak_file, sep='\t', nrows=1)
                self.logger.info(f"Header from file: {list(peaks_with_header.columns)}")
                
                # Re-read skipping the header line manually
                peaks = pd.read_csv(
                    peak_file, 
                    sep='\t', 
                    skiprows=1,
                    names=expected_columns,
                    low_memory=False
                )
            else:
                peaks.columns = expected_columns
            
            # Ensure numeric columns are properly typed
            numeric_cols = ['start', 'end', 'score', 'coverage', 'fold_change', 'pvalue', 'adjusted_pvalue', 'read_count']
            for col in numeric_cols:
                if col in peaks.columns:
                    peaks[col] = pd.to_numeric(peaks[col], errors='coerce')
            
            self.logger.info(f"Successfully loaded {len(peaks)} peaks")
            return peaks
            
        except Exception as e:
            self.logger.error(f"Error loading peaks from {peak_file}: {str(e)}")
            raise
    
    def validate_peak(self, chrom, start, end, fold_change=None, 
                     coverage=None, pvalue=None) -> CircularCandidate:
        """
        Comprehensive validation of a peak region
        Combines all evidence types to classify as eccDNA or artifact
        
        Parameters:
        -----------
        chrom : str
            Chromosome name
        start : int
            Peak start position
        end : int
            Peak end position
        fold_change : float, optional
            Fold change from coverage analysis
        coverage : float, optional
            Mean coverage in region
        pvalue : float, optional
            P-value from coverage analysis
            
        Returns:
        --------
        CircularCandidate
            Validated candidate with classification and evidence
        """
        
        self.logger.debug(f"Validating peak: {chrom}:{start}-{end}")
        
        # Initialize candidate
        candidate = CircularCandidate(
            chromosome=chrom,
            start=start,
            end=end,
            length=end - start,
            fold_change=fold_change or 0,
            coverage=coverage or 0,
            pvalue=pvalue or 1.0
        )
        
        try:
            # Step 1: Analyze read characteristics for artifact detection
            read_stats = self.analyze_read_characteristics(chrom, start, end)
            candidate.mean_mapq = read_stats['mean_mapq']
            candidate.soft_clip_rate = read_stats['soft_clip_rate']
            candidate.position_diversity = read_stats['position_diversity']
            
            # Step 2: Look for junction evidence supporting circular DNA
            junction_evidence = self.find_junction_evidence(chrom, start, end)
            candidate.junction_reads = junction_evidence['junction_reads']
            candidate.split_reads = junction_evidence['split_reads']
            candidate.discordant_pairs = junction_evidence['discordant_pairs']
            candidate.boundary_clips = junction_evidence['boundary_clips']
            
            # Step 3: Check coverage uniformity across the region
            coverage_stats = self.analyze_coverage_pattern(chrom, start, end)
            candidate.coverage_uniformity = coverage_stats['uniformity_score']
            
            # Step 4: Classify based on all evidence
            self.classify_candidate(candidate)
            
            # Step 5: Calculate confidence score
            self.calculate_confidence_score(candidate)
            
        except Exception as e:
            self.logger.warning(f"Error validating peak {chrom}:{start}-{end}: {str(e)}")
            candidate.classification = "ERROR"
            candidate.evidence_summary = f"Validation error: {str(e)}"
        
        return candidate
    
    def analyze_read_characteristics(self, chrom, start, end, max_reads=1000):
        """
        Analyze read characteristics to identify repetitive element artifacts
        
        This method examines read mapping quality, soft-clipping patterns,
        and position diversity to detect problematic regions.
        """
        
        stats = {
            'total_reads': 0,
            'unique_reads': set(),
            'mapq_scores': [],
            'soft_clipped': 0,
            'high_soft_clip': 0,  # >30% soft-clipped
            'extreme_soft_clip': 0,  # >50% soft-clipped
            'start_positions': [],
            'end_positions': []
        }
        
        try:
            for i, read in enumerate(self.bam.fetch(chrom, start, end)):
                if i >= max_reads:
                    break
                
                stats['total_reads'] += 1
                stats['unique_reads'].add(read.query_name)
                stats['mapq_scores'].append(read.mapping_quality)
                stats['start_positions'].append(read.reference_start)
                stats['end_positions'].append(read.reference_end)
                
                # Analyze CIGAR for soft-clipping patterns
                if read.cigartuples:
                    total_len = sum(l for op, l in read.cigartuples if op in [0, 1, 4, 5])
                    soft_clip_len = sum(l for op, l in read.cigartuples if op == 4)
                    
                    if soft_clip_len > 0:
                        stats['soft_clipped'] += 1
                        soft_clip_pct = soft_clip_len / total_len if total_len > 0 else 0
                        
                        if soft_clip_pct > 0.3:
                            stats['high_soft_clip'] += 1
                        if soft_clip_pct > 0.5:
                            stats['extreme_soft_clip'] += 1
        
        except Exception as e:
            self.logger.warning(f"Error analyzing reads in {chrom}:{start}-{end}: {str(e)}")
            return {
                'total_reads': 0, 'unique_read_rate': 0, 'mean_mapq': 0,
                'soft_clip_rate': 0, 'extreme_soft_clip_rate': 0, 'position_diversity': 0
            }
        
        # Calculate summary metrics
        return {
            'total_reads': stats['total_reads'],
            'unique_read_rate': len(stats['unique_reads']) / max(stats['total_reads'], 1),
            'mean_mapq': np.mean(stats['mapq_scores']) if stats['mapq_scores'] else 0,
            'soft_clip_rate': stats['high_soft_clip'] / max(stats['total_reads'], 1),
            'extreme_soft_clip_rate': stats['extreme_soft_clip'] / max(stats['total_reads'], 1),
            'position_diversity': len(set(stats['start_positions'])) / max(stats['total_reads'], 1)
        }
    
    def find_junction_evidence(self, chrom, start, end, window=500):
        """
        Look for junction reads that support circular DNA formation
        
        This method searches for:
        1. Split reads with supplementary alignments
        2. Discordant read pairs
        3. Soft-clipped reads at boundaries
        
        These patterns are characteristic of circular DNA junctions.
        """
        
        evidence = {
            'junction_reads': 0,
            'split_reads': 0,
            'discordant_pairs': 0,
            'boundary_clips': 0,
            'junction_positions': []
        }
        
        try:
            # Check reads near start boundary for junction evidence
            for read in self.bam.fetch(chrom, max(0, start - window), start + window):
                # Look for split reads with supplementary alignments (SA tag)
                if read.has_tag('SA'):
                    sa_tag = read.get_tag('SA')
                    sa_parts = sa_tag.split(',')
                    if len(sa_parts) >= 2:
                        try:
                            sa_chr = sa_parts[0]
                            sa_pos = int(sa_parts[1])
                            
                            # Check if supplementary alignment maps to the other end
                            if sa_chr == chrom and abs(sa_pos - end) < window:
                                evidence['split_reads'] += 1
                                evidence['junction_positions'].append((read.reference_start, sa_pos))
                        except (ValueError, IndexError):
                            pass
                
                # Check for soft-clipping at region boundaries
                if read.cigartuples and read.reference_start <= start + 50:
                    if read.cigartuples[-1][0] == 4:  # Right soft-clip
                        evidence['boundary_clips'] += 1
            
            # Check reads near end boundary
            try:
                ref_length = self.bam.get_reference_length(chrom)
                end_region_end = min(end + window, ref_length) if ref_length else end + window
            except:
                end_region_end = end + window
            
            for read in self.bam.fetch(chrom, end - window, end_region_end):
                if read.cigartuples and read.reference_end >= end - 50:
                    if read.cigartuples[0][0] == 4:  # Left soft-clip
                        evidence['boundary_clips'] += 1
                
                # Check for discordant pairs pointing to the other end
                if read.is_paired and not read.is_proper_pair:
                    if read.next_reference_name == chrom:
                        if abs(read.next_reference_start - start) < window:
                            evidence['discordant_pairs'] += 1
        
        except Exception as e:
            self.logger.warning(f"Error finding junction evidence in {chrom}:{start}-{end}: {str(e)}")
        
        # Calculate total junction support
        evidence['junction_reads'] = (evidence['split_reads'] + 
                                     evidence['discordant_pairs'] + 
                                     (evidence['boundary_clips'] // 2))
        
        return evidence
    
    def analyze_coverage_pattern(self, chrom, start, end, bin_size=100):
        """
        Analyze coverage uniformity across the region
        
        Uniform coverage suggests authentic eccDNA, while spiky or irregular
        coverage may indicate artifacts or repetitive elements.
        """
        
        num_bins = max(1, (end - start) // bin_size)
        bin_size = (end - start) // num_bins
        coverage_bins = []
        
        try:
            for i in range(num_bins):
                bin_start = start + i * bin_size
                bin_end = min(bin_start + bin_size, end)
                
                # Count reads in this bin
                read_count = self.bam.count(chrom, bin_start, bin_end)
                coverage_bins.append(read_count)
            
            if len(coverage_bins) > 1 and np.mean(coverage_bins) > 0:
                # Calculate coefficient of variation
                cv = np.std(coverage_bins) / np.mean(coverage_bins)
                uniformity_score = 1 / (1 + cv)
            else:
                uniformity_score = 0.5 if len(coverage_bins) == 1 else 0
        
        except Exception as e:
            self.logger.warning(f"Error analyzing coverage in {chrom}:{start}-{end}: {str(e)}")
            coverage_bins = []
            uniformity_score = 0
        
        return {
            'coverage_bins': coverage_bins,
            'mean_coverage': np.mean(coverage_bins) if coverage_bins else 0,
            'cv': cv if 'cv' in locals() else None,
            'uniformity_score': uniformity_score
        }
    
    def classify_candidate(self, candidate: CircularCandidate):
        """
        Classify candidate as eccDNA, artifact, or uncertain based on all evidence
        
        Uses CircONTrack-specific criteria for classification.
        """
        
        artifact_reasons = []
        eccDNA_evidence = []
        
        # Check for artifact signatures
        if candidate.fold_change > self.max_fold_change:
            artifact_reasons.append(f"Extreme fold change ({candidate.fold_change:.0f}x)")
        
        if candidate.mean_mapq < self.min_mapq:
            artifact_reasons.append(f"Low MAPQ ({candidate.mean_mapq:.1f})")
        
        if candidate.soft_clip_rate > self.max_soft_clip_rate:
            artifact_reasons.append(f"High soft-clipping ({candidate.soft_clip_rate:.1%})")
        
        if candidate.position_diversity < 0.1:
            artifact_reasons.append("Low read position diversity")
        
        if candidate.coverage_uniformity < 0.3:
            artifact_reasons.append("Non-uniform coverage")
        
        # Check for eccDNA supporting evidence
        if candidate.junction_reads >= self.min_junction_support:
            eccDNA_evidence.append(f"Junction support ({candidate.junction_reads} reads)")
        
        if candidate.split_reads > 0:
            eccDNA_evidence.append(f"Split reads ({candidate.split_reads})")
        
        if candidate.coverage_uniformity > 0.7:
            eccDNA_evidence.append("Uniform coverage")
        
        if 2 <= candidate.fold_change <= 50:
            eccDNA_evidence.append(f"Reasonable fold change ({candidate.fold_change:.1f}x)")
        
        if candidate.mean_mapq >= 30:
            eccDNA_evidence.append("High mapping quality")
        
        # Make classification decision
        if len(artifact_reasons) >= 2:
            candidate.classification = "ARTIFACT"
            candidate.is_artifact = True
            candidate.artifact_reasons = artifact_reasons
        elif len(eccDNA_evidence) >= 2:
            candidate.classification = "LIKELY_ECCDNA"
            candidate.is_artifact = False
        elif len(eccDNA_evidence) >= 1 and len(artifact_reasons) == 0:
            candidate.classification = "POSSIBLE_ECCDNA"
            candidate.is_artifact = False
        else:
            candidate.classification = "UNCERTAIN"
            candidate.is_artifact = False
        
        # Create evidence summary
        evidence_parts = []
        if eccDNA_evidence:
            evidence_parts.append("Supporting: " + "; ".join(eccDNA_evidence))
        if artifact_reasons:
            evidence_parts.append("Concerns: " + "; ".join(artifact_reasons))
        
        candidate.evidence_summary = " | ".join(evidence_parts) if evidence_parts else "No clear evidence"
    
    def calculate_confidence_score(self, candidate: CircularCandidate):
        """
        Calculate confidence score (0-100) based on all available evidence
        
        Higher scores indicate stronger evidence for authentic eccDNA.
        """
        
        score = 50  # Start with neutral score
        
        # Positive factors (eccDNA evidence)
        if candidate.junction_reads > 0:
            score += min(candidate.junction_reads * 5, 20)
        
        if candidate.split_reads > 0:
            score += min(candidate.split_reads * 5, 15)
        
        if candidate.coverage_uniformity > 0.7:
            score += 15
        elif candidate.coverage_uniformity > 0.5:
            score += 10
        
        if 2 <= candidate.fold_change <= 20:
            score += 10
        elif 20 < candidate.fold_change <= 50:
            score += 5
        
        if candidate.mean_mapq >= 30:
            score += 10
        elif candidate.mean_mapq >= 20:
            score += 5
        
        # Negative factors (artifact indicators)
        if candidate.fold_change > 100:
            score -= 30
        elif candidate.fold_change > 50:
            score -= 15
        
        if candidate.mean_mapq < 20:
            score -= 20
        elif candidate.mean_mapq < 30:
            score -= 10
        
        if candidate.soft_clip_rate > 0.5:
            score -= 25
        elif candidate.soft_clip_rate > 0.3:
            score -= 15
        
        if candidate.position_diversity < 0.1:
            score -= 20
        elif candidate.position_diversity < 0.3:
            score -= 10
        
        # Cap between 0 and 100
        candidate.confidence_score = max(0, min(100, score))
    
    def process_all_peaks(self, output_prefix="validated"):
        """
        Process and validate all peaks from the input peak file
        
        Parameters:
        -----------
        output_prefix : str
            Prefix for output files
            
        Returns:
        --------
        pd.DataFrame
            Results DataFrame with all validated peaks
        """
        
        if self.peaks is None:
            raise ValueError("No peaks loaded. Please provide a peak file.")
        
        self.logger.info(f"Processing {len(self.peaks)} peaks for validation...")
        
        validated_candidates = []
        
        for idx, peak in self.peaks.iterrows():
            if idx % 10 == 0:
                self.logger.info(f"Processing peak {idx+1}/{len(self.peaks)}")
            
            try:
                candidate = self.validate_peak(
                    peak['chr'], 
                    peak['start'], 
                    peak['end'],
                    peak['fold_change'],
                    peak['coverage'],
                    peak.get('adjusted_pvalue', peak.get('pvalue', 1.0))
                )
                
                validated_candidates.append(candidate)
                
            except Exception as e:
                self.logger.warning(f"Error processing peak {idx}: {str(e)}")
                # Create error candidate
                error_candidate = CircularCandidate(
                    chromosome=str(peak['chr']),
                    start=int(peak['start']),
                    end=int(peak['end']),
                    length=int(peak['end']) - int(peak['start']),
                    fold_change=float(peak['fold_change']),
                    coverage=float(peak['coverage']),
                    classification="ERROR",
                    evidence_summary=f"Processing error: {str(e)}"
                )
                validated_candidates.append(error_candidate)
        
        # Convert to DataFrame
        results_df = pd.DataFrame([{
            'chr': c.chromosome,
            'start': c.start,
            'end': c.end,
            'length': c.length,
            'fold_change': c.fold_change,
            'coverage': c.coverage,
            'pvalue': c.pvalue,
            'junction_reads': c.junction_reads,
            'split_reads': c.split_reads,
            'discordant_pairs': c.discordant_pairs,
            'boundary_clips': c.boundary_clips,
            'mean_mapq': c.mean_mapq,
            'soft_clip_rate': c.soft_clip_rate,
            'coverage_uniformity': c.coverage_uniformity,
            'confidence_score': c.confidence_score,
            'classification': c.classification,
            'is_artifact': c.is_artifact,
            'evidence': c.evidence_summary
        } for c in validated_candidates])
        
        # Save results
        self.save_results(results_df, output_prefix)
        
        return results_df
    
    def save_results(self, results_df, output_prefix):
        """Save validated results in multiple CircONTrack-compatible formats"""
        
        try:
            # Save full results with CircONTrack header
            full_output = f"{output_prefix}_all_peaks.tsv"
            with open(full_output, 'w') as f:
                f.write("# CircONTrack Peak Validation Results\n")
                f.write(f"# Total peaks processed: {len(results_df)}\n")
                f.write(f"# BAM file: {self.bam_file}\n")
                f.write(f"# Peak file: {self.peak_file}\n")
            
            results_df.to_csv(full_output, sep='\t', index=False, mode='a')
            self.logger.info(f"Saved all results to {full_output}")
            
            # Save likely eccDNA candidates
            eccdna_df = results_df[results_df['classification'].isin(['LIKELY_ECCDNA', 'POSSIBLE_ECCDNA'])]
            if len(eccdna_df) > 0:
                eccdna_output = f"{output_prefix}_eccdna.tsv"
                eccdna_df.to_csv(eccdna_output, sep='\t', index=False)
                self.logger.info(f"Saved {len(eccdna_df)} eccDNA candidates to {eccdna_output}")
                
                # BED format for genome browsers
                bed_output = f"{output_prefix}_eccdna.bed"
                bed_df = eccdna_df[['chr', 'start', 'end']].copy()
                bed_df['name'] = eccdna_df['classification']
                bed_df['score'] = eccdna_df['confidence_score'].astype(int)
                bed_df['strand'] = '.'
                bed_df.to_csv(bed_output, sep='\t', index=False, header=False)
                self.logger.info(f"Saved BED format to {bed_output}")
            
            # Save artifacts for inspection
            artifact_df = results_df[results_df['is_artifact']]
            if len(artifact_df) > 0:
                artifact_output = f"{output_prefix}_artifacts.tsv"
                artifact_df.to_csv(artifact_output, sep='\t', index=False)
                self.logger.info(f"Identified {len(artifact_df)} artifacts saved to {artifact_output}")
            
            # Print summary to console
            self.print_summary(results_df)
            
        except Exception as e:
            self.logger.error(f"Error saving results: {str(e)}")
            raise
    
    def print_summary(self, results_df):
        """Print comprehensive analysis summary"""
        
        print("\n" + "="*70)
        print("CIRCONTRACK PEAK VALIDATION SUMMARY")
        print(f"Version: {__version__}")
        print("="*70)
        
        print(f"\nTotal peaks analyzed: {len(results_df)}")
        
        # Classification breakdown
        class_counts = results_df['classification'].value_counts()
        print("\nClassification Results:")
        for class_name, count in class_counts.items():
            pct = count / len(results_df) * 100
            print(f"  {class_name}: {count} ({pct:.1f}%)")
        
        # High-confidence eccDNA
        high_conf = results_df[
            (results_df['classification'] == 'LIKELY_ECCDNA') & 
            (results_df['confidence_score'] >= 70)
        ]
        
        if len(high_conf) > 0:
            print(f"\nHigh-confidence eccDNA (score ≥70): {len(high_conf)}")
            print("\nTop 5 eccDNA candidates:")
            for _, row in high_conf.nlargest(5, 'confidence_score').iterrows():
                print(f"  {row['chr']}:{row['start']}-{row['end']}")
                print(f"    Confidence: {row['confidence_score']:.0f}/100")
                print(f"    Fold change: {row['fold_change']:.1f}x")
                print(f"    Junction support: {row['junction_reads']}")
                print(f"    Evidence: {row['evidence']}")
                print()
        
        # Artifact analysis
        artifacts = results_df[results_df['is_artifact']]
        if len(artifacts) > 0:
            print(f"\nArtifacts identified: {len(artifacts)}")
            
            # Check for extreme cases
            extreme = artifacts[artifacts['fold_change'] > 500]
            if len(extreme) > 0:
                print(f"  Extreme artifacts (>500x fold change): {len(extreme)}")
        
        # Junction evidence summary
        with_junctions = results_df[results_df['junction_reads'] > 0]
        if len(with_junctions) > 0:
            print(f"\nPeaks with junction evidence: {len(with_junctions)}")
            print(f"  Mean junction reads: {with_junctions['junction_reads'].mean():.1f}")
            print(f"  Max junction reads: {with_junctions['junction_reads'].max()}")
        
        print("\n" + "="*70)


def main():
    """Main entry point for CircONTrack peak validation"""
    
    parser = argparse.ArgumentParser(
        description="CircONTrack Peak Validation - Validate peaks with junction detection and artifact filtering",
        prog="circontrack-validate"
    )
    
    parser.add_argument('peak_file', 
                       help='Peak file from CircONTrack coverage peak caller')
    parser.add_argument('bam_file', 
                       help='BAM file used for peak calling')
    parser.add_argument('-r', '--reference', 
                       help='Reference FASTA file (optional)')
    parser.add_argument('-o', '--output-prefix', default='validated',
                       help='Output file prefix (default: validated)')
    parser.add_argument('--min-junction', type=int, default=2,
                       help='Minimum junction reads for eccDNA classification (default: 2)')
    parser.add_argument('--max-fold', type=float, default=100,
                       help='Maximum fold change before artifact flag (default: 100)')
    parser.add_argument('--min-mapq', type=float, default=20,
                       help='Minimum mean MAPQ score (default: 20)')
    parser.add_argument('--max-softclip', type=float, default=0.3,
                       help='Maximum soft-clipping rate (default: 0.3)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(verbose=args.verbose)
    
    try:
        logger.info("Starting CircONTrack peak validation...")
        
        # Initialize validation pipeline
        validator = CircONTrackValidator(
            bam_file=args.bam_file,
            peak_file=args.peak_file,
            reference_file=args.reference,
            min_junction_support=args.min_junction,
            max_fold_change=args.max_fold,
            min_mapq=args.min_mapq,
            max_soft_clip_rate=args.max_softclip,
            logger=logger
        )
        
        # Process all peaks
        logger.info("Running validation pipeline...")
        results = validator.process_all_peaks(args.output_prefix)
        
        # Summary statistics
        total_peaks = len(results)
        eccdna_candidates = len(results[results['classification'].isin(['LIKELY_ECCDNA', 'POSSIBLE_ECCDNA'])])
        artifacts = len(results[results['is_artifact']])
        
        logger.info(f"Validation complete!")
        logger.info(f"Results: {eccdna_candidates} eccDNA candidates, {artifacts} artifacts from {total_peaks} peaks")
        logger.info(f"Output files saved with prefix: {args.output_prefix}")
        
    except Exception as e:
        logger.error(f"Validation failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()