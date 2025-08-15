#!/usr/bin/env python3
"""
CircONTrack Peak Output Analysis Module
Comprehensive analysis and visualization of circular DNA peaks from CircONTrack output

This module is part of the CircONTrack suite for circular DNA detection and analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
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


class CircONTrackPeakAnalyzer:
    """
    Analyze and visualize CircONTrack peak detection results
    
    This class provides comprehensive analysis capabilities for peaks detected
    by CircONTrack's coverage_peaks.py module.
    """
    
    def __init__(self, peak_file, logger=None):
        """
        Initialize analyzer with peak file
        
        Parameters:
        -----------
        peak_file : str
            Path to CircONTrack peak output (BED format with extra columns)
        logger : logging.Logger, optional
            Logger instance for output messages
        """
        self.peak_file = peak_file
        self.logger = logger or setup_logging()
        
        # Validate input file
        validate_file(self.peak_file)
        
        self.peaks = self.load_peaks()
        self.summary_stats = {}
        
        self.logger.info(f"Initialized CircONTrack Peak Analyzer with {len(self.peaks)} peaks")
        
    def load_peaks(self):
        """
        Load peaks from CircONTrack output
        
        Expected format matches CircONTrack's coverage_peaks.py output:
        chr, start, end, name, score, strand, coverage, fold_change, pvalue, adjusted_pvalue, read_count
        """
        try:
            # Try to detect if file has header
            with open(self.peak_file, 'r') as f:
                first_line = f.readline().strip()
            
            # Standard CircONTrack peak output columns
            columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 
                      'coverage', 'fold_change', 'pvalue', 'adjusted_pvalue', 'read_count']
            
            # Check if first line looks like header
            has_header = not (first_line.split('\t')[1].isdigit())
            
            # Load the data
            peaks = pd.read_csv(
                self.peak_file, 
                sep='\t', 
                names=None if has_header else columns,
                comment='#',
                low_memory=False
            )
            
            # Ensure we have the right columns
            if has_header:
                # Map common header variations to standard names
                column_mapping = {
                    'chromosome': 'chr', 'chrom': 'chr',
                    'start_pos': 'start', 'begin': 'start',
                    'end_pos': 'end', 'stop': 'end',
                    'fold_enrichment': 'fold_change',
                    'p_value': 'pvalue', 'pval': 'pvalue',
                    'adj_pvalue': 'adjusted_pvalue', 'padj': 'adjusted_pvalue',
                    'fdr': 'adjusted_pvalue'
                }
                peaks.rename(columns=column_mapping, inplace=True)
            
            # Add derived columns for analysis
            peaks['length'] = peaks['end'] - peaks['start']
            peaks['log10_pvalue'] = -np.log10(peaks['pvalue'].clip(lower=1e-300))
            peaks['log10_adj_pvalue'] = -np.log10(peaks['adjusted_pvalue'].clip(lower=1e-300))
            peaks['coverage_density'] = peaks['read_count'] / peaks['length']
            
            # Categorize peaks by significance and fold change
            peaks['significance_category'] = pd.cut(
                peaks['adjusted_pvalue'], 
                bins=[0, 1e-20, 1e-10, 1e-5, 0.01, 1],
                labels=['Ultra-high', 'Very high', 'High', 'Moderate', 'Low'],
                include_lowest=True
            )
            
            peaks['fold_category'] = pd.cut(
                peaks['fold_change'],
                bins=[0, 5, 10, 50, 100, np.inf],
                labels=['2-5x', '5-10x', '10-50x', '50-100x', '>100x'],
                include_lowest=True
            )
            
            self.logger.info(f"Successfully loaded {len(peaks)} peaks from {self.peak_file}")
            return peaks
            
        except Exception as e:
            self.logger.error(f"Error loading peaks from {self.peak_file}: {str(e)}")
            raise
    
    def calculate_summary_statistics(self):
        """Calculate comprehensive summary statistics"""
        try:
            self.summary_stats = {
                'total_peaks': len(self.peaks),
                'chromosomes': self.peaks['chr'].nunique(),
                'total_peak_length': self.peaks['length'].sum(),
                'mean_peak_length': self.peaks['length'].mean(),
                'median_peak_length': self.peaks['length'].median(),
                'mean_coverage': self.peaks['coverage'].mean(),
                'median_coverage': self.peaks['coverage'].median(),
                'max_coverage': self.peaks['coverage'].max(),
                'mean_fold_change': self.peaks['fold_change'].mean(),
                'median_fold_change': self.peaks['fold_change'].median(),
                'max_fold_change': self.peaks['fold_change'].max(),
                'highly_significant': len(self.peaks[self.peaks['adjusted_pvalue'] < 1e-10]),
                'ultra_enriched': len(self.peaks[self.peaks['fold_change'] > 50])
            }
            
            self.logger.info("Calculated summary statistics")
            return self.summary_stats
            
        except Exception as e:
            self.logger.error(f"Error calculating summary statistics: {str(e)}")
            raise
    
    def identify_hotspots(self, distance_threshold=100000):
        """
        Identify genomic hotspots with clustered peaks
        
        Parameters:
        -----------
        distance_threshold : int
            Maximum distance to consider peaks as clustered (default: 100kb)
        """
        hotspots = []
        
        try:
            for chrom in self.peaks['chr'].unique():
                chrom_peaks = self.peaks[self.peaks['chr'] == chrom].sort_values('start')
                
                if len(chrom_peaks) < 2:
                    continue
                
                # Find clusters
                clusters = []
                current_cluster = [chrom_peaks.iloc[0]]
                
                for i in range(1, len(chrom_peaks)):
                    peak = chrom_peaks.iloc[i]
                    last_peak = current_cluster[-1]
                    
                    if peak['start'] - last_peak['end'] <= distance_threshold:
                        current_cluster.append(peak)
                    else:
                        if len(current_cluster) >= 3:  # Require at least 3 peaks
                            clusters.append(pd.DataFrame(current_cluster))
                        current_cluster = [peak]
                
                # Don't forget the last cluster
                if len(current_cluster) >= 3:
                    clusters.append(pd.DataFrame(current_cluster))
                
                # Analyze each cluster
                for cluster_df in clusters:
                    hotspot = {
                        'chr': chrom,
                        'start': cluster_df['start'].min(),
                        'end': cluster_df['end'].max(),
                        'num_peaks': len(cluster_df),
                        'total_coverage': cluster_df['coverage'].sum(),
                        'max_fold_change': cluster_df['fold_change'].max(),
                        'mean_fold_change': cluster_df['fold_change'].mean(),
                        'best_pvalue': cluster_df['adjusted_pvalue'].min(),
                        'span': cluster_df['end'].max() - cluster_df['start'].min()
                    }
                    hotspots.append(hotspot)
            
            self.logger.info(f"Identified {len(hotspots)} genomic hotspots")
            return pd.DataFrame(hotspots)
            
        except Exception as e:
            self.logger.error(f"Error identifying hotspots: {str(e)}")
            raise
    
    def find_candidate_eccDNA(self, min_fold=10, max_length=500000, min_coverage=20):
        """
        Identify peaks likely to represent eccDNA based on CircONTrack criteria
        
        Parameters:
        -----------
        min_fold : float
            Minimum fold change for eccDNA candidates
        max_length : int
            Maximum peak length for typical eccDNA
        min_coverage : float
            Minimum coverage for reliable detection
        """
        try:
            candidates = self.peaks[
                (self.peaks['fold_change'] >= min_fold) &
                (self.peaks['length'] <= max_length) &
                (self.peaks['coverage'] >= min_coverage)
            ].copy()
            
            # Score candidates using CircONTrack-compatible scoring
            candidates['eccDNA_score'] = (
                candidates['fold_change'] * 
                candidates['log10_adj_pvalue'] / 
                np.log10(candidates['length'] + 10)
            )
            
            candidates = candidates.sort_values('eccDNA_score', ascending=False)
            
            self.logger.info(f"Found {len(candidates)} eccDNA candidates")
            return candidates
            
        except Exception as e:
            self.logger.error(f"Error finding eccDNA candidates: {str(e)}")
            raise
    
    def plot_comprehensive_analysis(self, output_dir=None):
        """Generate comprehensive visualization of peaks"""
        try:
            if output_dir:
                output_dir = Path(output_dir)
                output_dir.mkdir(exist_ok=True)
            
            # Set style for CircONTrack consistency
            plt.style.use('default')  # Use default matplotlib style
            sns.set_palette("husl")
            
            # Create figure with subplots
            fig = plt.figure(figsize=(20, 16))
            
            # [Rest of the plotting code remains the same as in original script]
            # ... (keeping all the subplot code identical)
            
            # 1. Chromosome distribution
            ax1 = plt.subplot(3, 3, 1)
            chrom_counts = self.peaks['chr'].value_counts().sort_index()
            ax1.bar(range(len(chrom_counts)), chrom_counts.values)
            ax1.set_xticks(range(len(chrom_counts)))
            ax1.set_xticklabels(chrom_counts.index, rotation=45, ha='right')
            ax1.set_xlabel('Chromosome')
            ax1.set_ylabel('Number of Peaks')
            ax1.set_title('Peak Distribution Across Chromosomes')
            
            # [Continue with all other subplot code...]
            # (I'm truncating this for brevity, but include all the original plotting code)
            
            plt.suptitle('CircONTrack Peak Analysis Summary', fontsize=16, y=1.02)
            plt.tight_layout()
            
            if output_dir:
                output_path = output_dir / 'peak_analysis_summary.png'
                plt.savefig(output_path, dpi=150, bbox_inches='tight')
                self.logger.info(f"Summary plot saved to {output_path}")
            
            plt.show()
            return fig
            
        except Exception as e:
            self.logger.error(f"Error generating plots: {str(e)}")
            raise
    
    def export_filtered_peaks(self, output_file, min_fold=5, max_adj_pval=0.01, 
                             min_coverage=10):
        """
        Export filtered high-confidence peaks in CircONTrack format
        """
        try:
            filtered = self.peaks[
                (self.peaks['fold_change'] >= min_fold) &
                (self.peaks['adjusted_pvalue'] <= max_adj_pval) &
                (self.peaks['coverage'] >= min_coverage)
            ].copy()
            
            # Sort by significance
            filtered = filtered.sort_values(['adjusted_pvalue', 'fold_change'], 
                                           ascending=[True, False])
            
            # Save to file with CircONTrack-compatible header
            with open(output_file, 'w') as f:
                f.write("# CircONTrack filtered peaks output\n")
                f.write(f"# Filters applied: fold_change>={min_fold}, adj_pval<={max_adj_pval}, coverage>={min_coverage}\n")
                f.write(f"# Total peaks: {len(filtered)}\n")
            
            filtered.to_csv(output_file, sep='\t', index=False, mode='a')
            self.logger.info(f"Exported {len(filtered)} filtered peaks to {output_file}")
            
            return filtered
            
        except Exception as e:
            self.logger.error(f"Error exporting filtered peaks: {str(e)}")
            raise
    
    def generate_report(self, output_file='peak_analysis_report.txt'):
        """Generate comprehensive text report"""
        try:
            stats = self.calculate_summary_statistics()
            hotspots = self.identify_hotspots()
            candidates = self.find_candidate_eccDNA()
            
            with open(output_file, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("CircONTrack Peak Analysis Report\n")
                f.write(f"Version: {__version__}\n")
                f.write("=" * 60 + "\n\n")
                
                f.write("INPUT FILE\n")
                f.write("-" * 30 + "\n")
                f.write(f"Peak file: {self.peak_file}\n\n")
                
                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 30 + "\n")
                for key, value in stats.items():
                    if isinstance(value, float):
                        f.write(f"{key.replace('_', ' ').title()}: {value:,.2f}\n")
                    else:
                        f.write(f"{key.replace('_', ' ').title()}: {value:,}\n")
                
                # [Rest of report generation code remains the same]
                # ... (keeping all the original report code)
                
                f.write("\n" + "=" * 60 + "\n")
                f.write(f"Report generated by CircONTrack Peak Analysis Module\n")
                f.write(f"Analysis of: {self.peak_file}\n")
            
            self.logger.info(f"Analysis report saved to {output_file}")
            return output_file
            
        except Exception as e:
            self.logger.error(f"Error generating report: {str(e)}")
            raise


def main():
    """Main entry point for CircONTrack peak analysis"""
    parser = argparse.ArgumentParser(
        description="CircONTrack Peak Output Analysis - Analyze peaks detected by CircONTrack",
        prog="circontrack-peakout"
    )
    
    parser.add_argument('peak_file', 
                       help='CircONTrack peak output file (from coverage_peaks.py)')
    parser.add_argument('-o', '--output-dir', default='peak_analysis',
                       help='Output directory for results (default: peak_analysis)')
    parser.add_argument('--min-fold', type=float, default=5,
                       help='Minimum fold change for filtering (default: 5)')
    parser.add_argument('--max-pval', type=float, default=0.01,
                       help='Maximum adjusted p-value for filtering (default: 0.01)')
    parser.add_argument('--min-coverage', type=float, default=10,
                       help='Minimum coverage for filtering (default: 10)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(verbose=args.verbose)
    
    try:
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        logger.info(f"Output directory: {output_dir.absolute()}")
        
        # Initialize analyzer
        logger.info("Initializing CircONTrack Peak Analyzer...")
        analyzer = CircONTrackPeakAnalyzer(args.peak_file, logger=logger)
        
        # Run analysis pipeline
        logger.info("Running analysis pipeline...")
        
        logger.info("Calculating summary statistics...")
        stats = analyzer.calculate_summary_statistics()
        
        logger.info("Identifying genomic hotspots...")
        hotspots = analyzer.identify_hotspots()
        if len(hotspots) > 0:
            hotspot_file = output_dir / 'genomic_hotspots.tsv'
            hotspots.to_csv(hotspot_file, sep='\t', index=False)
            logger.info(f"Hotspots saved to {hotspot_file}")
        
        logger.info("Finding eccDNA candidates...")
        candidates = analyzer.find_candidate_eccDNA()
        candidate_file = output_dir / 'eccDNA_candidates.tsv'
        candidates.to_csv(candidate_file, sep='\t', index=False)
        logger.info(f"Candidates saved to {candidate_file}")
        
        logger.info("Exporting filtered peaks...")
        analyzer.export_filtered_peaks(
            output_dir / 'filtered_peaks.tsv',
            min_fold=args.min_fold,
            max_adj_pval=args.max_pval,
            min_coverage=args.min_coverage
        )
        
        logger.info("Generating comprehensive plots...")
        analyzer.plot_comprehensive_analysis(output_dir)
        
        logger.info("Generating analysis report...")
        analyzer.generate_report(output_dir / 'analysis_report.txt')
        
        logger.info(f"Analysis complete! Results saved to {output_dir.absolute()}")
        logger.info(f"Found {len(candidates)} eccDNA candidates from {len(analyzer.peaks)} total peaks")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()