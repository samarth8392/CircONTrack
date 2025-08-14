#!/usr/bin/env python3
"""
CircONTrack Coverage Peaks - Statistical coverage peak detection using negative binomial distribution
Identifies statistically significant coverage enrichments in ONT sequencing data
"""

import numpy as np
import pysam
from scipy import stats
from scipy.signal import find_peaks
from statsmodels.stats.multitest import multipletests
import argparse
import sys
import matplotlib
matplotlib.use('Agg') 
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Optional imports for visualization
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False

class CoveragePeakCaller:
    """
    Statistical coverage peak detection using negative binomial distribution
    Designed for ONT long-read data with overdispersed coverage patterns
    """
    
    def __init__(self, bam_file, reference_file=None, window_size=1000, step_size=None):
        """
        Initialize coverage peak caller
        
        Parameters:
        -----------
        bam_file : str
            Path to indexed BAM file
        reference_file : str
            Path to reference FASTA (optional, for getting chromosome lengths)
        window_size : int
            Size of sliding window in bp (default: 1000)
        step_size : int
            Step size for sliding window (default: window_size/2)
        """
        self.bam_file = bam_file
        self.reference_file = reference_file
        self.window_size = window_size
        self.step_size = step_size or window_size // 2
        
        # Open BAM file
        self.bam = pysam.AlignmentFile(bam_file, "rb")
        
        # Get chromosome lengths
        self.chrom_lengths = {}
        if reference_file:
            ref = pysam.FastaFile(reference_file)
            self.chrom_lengths = {name: ref.get_reference_length(name) 
                                 for name in ref.references}
            ref.close()
        else:
            # Get from BAM header
            for i, name in enumerate(self.bam.references):
                self.chrom_lengths[name] = self.bam.lengths[i]
        
        print(f"Initialized peak caller with {len(self.chrom_lengths)} references")
        print(f"Window size: {window_size} bp, Step size: {step_size} bp")
    
    def calculate_coverage_windows(self, chromosome, start=None, end=None):
        """
        Calculate read coverage in sliding windows
        
        Returns:
        --------
        windows : list of dict
            Each window contains: chr, start, end, coverage, read_count
        """
        if chromosome not in self.chrom_lengths:
            raise ValueError(f"Chromosome {chromosome} not found in reference")
        
        chrom_length = self.chrom_lengths[chromosome]
        start = start or 0
        end = end or chrom_length
        
        windows = []
        
        # Slide windows across region
        for window_start in range(start, end, self.step_size):
            window_end = min(window_start + self.window_size, end)
            
            # Count reads in window
            read_count = 0
            coverage_sum = 0
            
            for read in self.bam.fetch(chromosome, window_start, window_end):
                if read.is_unmapped or read.is_secondary:
                    continue
                
                # Calculate overlap with window
                read_start = max(read.reference_start, window_start)
                read_end = min(read.reference_end, window_end)
                overlap = read_end - read_start
                
                if overlap > 0:
                    read_count += 1
                    coverage_sum += overlap
            
            # Calculate average coverage
            window_length = window_end - window_start
            avg_coverage = coverage_sum / window_length if window_length > 0 else 0
            
            windows.append({
                'chr': chromosome,
                'start': window_start,
                'end': window_end,
                'coverage': avg_coverage,
                'read_count': read_count,
                'length': window_length
            })
        
        return windows
    
    def fit_negative_binomial(self, coverage_values):
        """
        Fit negative binomial distribution to coverage data
        
        The negative binomial is appropriate for overdispersed count data
        (variance > mean), which is common in sequencing coverage
        
        Returns:
        --------
        params : dict
            Contains 'n' (number of failures) and 'p' (probability)
            Also 'mean' and 'dispersion' for interpretation
        """
        coverage_values = np.array(coverage_values)
        
        # Remove zeros for better fit (can handle separately)
        non_zero = coverage_values[coverage_values > 0]
        
        if len(non_zero) < 10:
            raise ValueError("Insufficient non-zero coverage values for fitting")
        
        # Method of moments estimation
        mean_cov = np.mean(non_zero)
        var_cov = np.var(non_zero)
        
        # Check for overdispersion
        if var_cov <= mean_cov:
            # Use Poisson if not overdispersed
            print("Warning: Data not overdispersed, using Poisson approximation")
            return {
                'distribution': 'poisson',
                'mean': mean_cov,
                'lambda': mean_cov,
                'dispersion': 1.0
            }
        
        # Negative binomial parameters
        # Mean = n(1-p)/p, Variance = n(1-p)/p^2
        p = mean_cov / var_cov
        n = mean_cov * p / (1 - p)
        
        # Calculate dispersion parameter (alpha)
        # Dispersion = 1/n (smaller n = more dispersion)
        dispersion = 1 / n if n > 0 else float('inf')
        
        return {
            'distribution': 'nbinom',
            'n': n,
            'p': p,
            'mean': mean_cov,
            'variance': var_cov,
            'dispersion': dispersion
        }
    
    def calculate_pvalues(self, windows, dist_params, use_median=False):
        """
        Calculate p-values for each window using fitted distribution
        
        Parameters:
        -----------
        windows : list of dict
            Coverage windows
        dist_params : dict
            Parameters from fit_negative_binomial
        use_median : bool
            Use median instead of mean for baseline
        
        Returns:
        --------
        pvalues : array
            P-value for each window (one-tailed test for enrichment)
        """
        coverage_values = np.array([w['coverage'] for w in windows])
        
        if use_median:
            baseline = np.median(coverage_values[coverage_values > 0])
        else:
            baseline = dist_params['mean']
        
        pvalues = []
        
        for window in windows:
            coverage = window['coverage']
            
            if dist_params['distribution'] == 'nbinom':
                # Negative binomial test
                # P(X >= observed | H0: coverage from background distribution)
                # Use survival function (1 - CDF)
                pval = stats.nbinom.sf(coverage - 1, 
                                      dist_params['n'], 
                                      dist_params['p'])
            else:
                # Poisson test
                pval = stats.poisson.sf(coverage - 1, 
                                       dist_params['lambda'])
            
            pvalues.append(pval)
        
        return np.array(pvalues)
    
    def call_peaks(self, chromosome, fdr_threshold=0.05, min_fold_change=2.0,
                   min_coverage=5, merge_distance=1000, plot=False):
        """
        Call statistically significant coverage peaks
        
        Parameters:
        -----------
        chromosome : str
            Chromosome to analyze
        fdr_threshold : float
            False discovery rate threshold (default: 0.05)
        min_fold_change : float
            Minimum fold change over background (default: 2.0)
        min_coverage : float
            Minimum average coverage in peak (default: 5)
        merge_distance : int
            Merge peaks within this distance (default: 1000 bp)
        plot : bool
            Generate coverage plot with peaks
        
        Returns:
        --------
        peaks : list of dict
            Significant peaks with statistics
        """
        print(f"\nAnalyzing {chromosome}...")
        
        # Calculate coverage windows
        windows = self.calculate_coverage_windows(chromosome)
        
        if len(windows) < 10:
            print(f"  Insufficient windows for {chromosome}")
            return []
        
        coverage_values = np.array([w['coverage'] for w in windows])
        
        # Fit distribution
        try:
            dist_params = self.fit_negative_binomial(coverage_values)
            print(f"  Fitted {dist_params['distribution']} distribution:")
            print(f"    Mean coverage: {dist_params['mean']:.2f}")
            print(f"    Variance: {dist_params.get('variance', 0):.2f}")
            print(f"    Dispersion: {dist_params.get('dispersion', 0):.4f}")
        except Exception as e:
            print(f"  Error fitting distribution: {e}")
            return []
        
        # Calculate p-values
        pvalues = self.calculate_pvalues(windows, dist_params)
        
        # Multiple testing correction
        if len(pvalues) > 0:
            # Benjamini-Hochberg FDR correction
            rejected, adjusted_pvals, _, _ = multipletests(
                pvalues, 
                alpha=fdr_threshold, 
                method='fdr_bh'
            )
            
            # Alternative: Bonferroni correction (more conservative)
            # rejected_bonf, adjusted_pvals_bonf, _, _ = multipletests(
            #     pvalues, 
            #     alpha=fdr_threshold, 
            #     method='bonferroni'
            # )
        else:
            return []
        
        # Filter for significant peaks
        peaks = []
        background = dist_params['mean']
        
        for i, (window, pval, adj_pval, is_sig) in enumerate(
            zip(windows, pvalues, adjusted_pvals, rejected)):
            
            fold_change = window['coverage'] / background if background > 0 else 0
            
            # Apply filters
            if (is_sig and 
                fold_change >= min_fold_change and 
                window['coverage'] >= min_coverage):
                
                peaks.append({
                    'chr': window['chr'],
                    'start': window['start'],
                    'end': window['end'],
                    'coverage': window['coverage'],
                    'read_count': window['read_count'],
                    'pvalue': pval,
                    'adjusted_pvalue': adj_pval,
                    'fold_change': fold_change,
                    'window_index': i
                })
        
        print(f"  Windows tested: {len(windows)}")
        print(f"  Significant windows (FDR < {fdr_threshold}): {sum(rejected)}")
        print(f"  Peaks after filtering: {len(peaks)}")
        
        # Merge nearby peaks
        if peaks and merge_distance > 0:
            peaks = self.merge_peaks(peaks, merge_distance)
            print(f"  Peaks after merging: {len(peaks)}")
        
        # Plot if requested
        if plot and HAS_PLOTTING and peaks:
            self.plot_coverage_peaks(chromosome, windows, peaks, dist_params)
        
        return peaks
    
    def merge_peaks(self, peaks, max_distance):
        """
        Merge peaks within specified distance
        
        Parameters:
        -----------
        peaks : list of dict
            Peaks to merge
        max_distance : int
            Maximum distance between peaks to merge
        
        Returns:
        --------
        merged : list of dict
            Merged peaks with combined statistics
        """
        if not peaks:
            return []
        
        # Sort by position
        peaks = sorted(peaks, key=lambda x: (x['chr'], x['start']))
        
        merged = []
        current = peaks[0].copy()
        
        for peak in peaks[1:]:
            if (peak['chr'] == current['chr'] and 
                peak['start'] - current['end'] <= max_distance):
                
                # Merge peaks
                current['end'] = peak['end']
                current['coverage'] = max(current['coverage'], peak['coverage'])
                current['read_count'] += peak['read_count']
                current['pvalue'] = min(current['pvalue'], peak['pvalue'])
                current['adjusted_pvalue'] = min(current['adjusted_pvalue'], 
                                                peak['adjusted_pvalue'])
                current['fold_change'] = max(current['fold_change'], 
                                            peak['fold_change'])
            else:
                # Save current and start new
                merged.append(current)
                current = peak.copy()
        
        merged.append(current)
        
        # Update peak lengths
        for peak in merged:
            peak['length'] = peak['end'] - peak['start']
        
        return merged
    
    def call_peaks_genome_wide(self, chromosomes=None, **kwargs):
        """
        Call peaks across all chromosomes
        
        Parameters:
        -----------
        chromosomes : list
            Specific chromosomes to analyze (default: all)
        **kwargs : 
            Arguments passed to call_peaks()
        
        Returns:
        --------
        all_peaks : dict
            Peaks organized by chromosome
        """
        if chromosomes is None:
            chromosomes = list(self.chrom_lengths.keys())
        
        all_peaks = {}
        total_peaks = 0
        
        for chrom in chromosomes:
            peaks = self.call_peaks(chrom, **kwargs)
            if peaks:
                all_peaks[chrom] = peaks
                total_peaks += len(peaks)
        
        print(f"\nTotal peaks called: {total_peaks}")
        return all_peaks
    
    def write_bed(self, peaks, output_file, name_prefix="peak"):
        """
        Write peaks to BED format
        
        Parameters:
        -----------
        peaks : dict or list
            Peaks to write (dict by chromosome or flat list)
        output_file : str
            Output BED file path
        name_prefix : str
            Prefix for peak names
        """
        with open(output_file, 'w') as f:
            # Header
            f.write("#chr\tstart\tend\tname\tscore\tstrand\tcoverage\tfold_change\t")
            f.write("pvalue\tadjusted_pvalue\tread_count\n")
            
            peak_id = 1
            
            # Handle both dict and list formats
            if isinstance(peaks, dict):
                # Dictionary organized by chromosome
                for chrom in sorted(peaks.keys()):
                    for peak in peaks[chrom]:
                        self._write_peak_line(f, peak, peak_id, name_prefix)
                        peak_id += 1
            else:
                # Flat list of peaks
                for peak in peaks:
                    self._write_peak_line(f, peak, peak_id, name_prefix)
                    peak_id += 1
        
        print(f"Peaks written to: {output_file}")
    
    def _write_peak_line(self, file_handle, peak, peak_id, name_prefix):
        """Write a single peak line to BED file"""
        # Convert p-value to phred-like score (0-1000)
        score = min(int(-10 * np.log10(peak['pvalue'] + 1e-300)), 1000)
        
        file_handle.write(f"{peak['chr']}\t{peak['start']}\t{peak['end']}\t")
        file_handle.write(f"{name_prefix}_{peak_id}\t{score}\t.\t")
        file_handle.write(f"{peak['coverage']:.2f}\t{peak['fold_change']:.2f}\t")
        file_handle.write(f"{peak['pvalue']:.2e}\t{peak['adjusted_pvalue']:.2e}\t")
        file_handle.write(f"{peak['read_count']}\n")
    
    def plot_coverage_peaks(self, chromosome, windows, peaks, dist_params, 
                           output_file=None):
        """
        Plot coverage profile with called peaks
        
        Parameters:
        -----------
        chromosome : str
            Chromosome name
        windows : list
            Coverage windows
        peaks : list
            Called peaks
        dist_params : dict
            Distribution parameters
        output_file : str
            Save plot to file (optional)
        """
        if not HAS_PLOTTING:
            print("Matplotlib not available for plotting")
            return
        
        # Extract data
        positions = np.array([w['start'] for w in windows])
        coverage = np.array([w['coverage'] for w in windows])
        
        # Create figure
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10), 
                                            height_ratios=[3, 1, 1])
        
        # Main coverage plot
        ax1.plot(positions, coverage, 'b-', alpha=0.6, linewidth=0.5, 
                label='Coverage')
        ax1.axhline(y=dist_params['mean'], color='gray', linestyle='--', 
                   alpha=0.5, label=f"Mean: {dist_params['mean']:.1f}")
        
        # Highlight peaks
        for peak in peaks:
            ax1.axvspan(peak['start'], peak['end'], alpha=0.3, color='red')
        
        ax1.set_ylabel('Coverage (X)')
        ax1.set_title(f'Coverage Profile - {chromosome}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # P-value track
        pvalues = []
        for w in windows:
            # Find if window is in a peak
            in_peak = False
            for peak in peaks:
                if w['start'] >= peak['start'] and w['end'] <= peak['end']:
                    pval = peak['pvalue']
                    in_peak = True
                    break
            if not in_peak:
                pval = 1.0
            pvalues.append(pval)
        
        pvalues = np.array(pvalues)
        log_pvals = -np.log10(pvalues + 1e-300)
        
        ax2.plot(positions, log_pvals, 'g-', alpha=0.6, linewidth=0.5)
        ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', 
                   alpha=0.5, label='FDR = 0.05')
        ax2.set_ylabel('-log10(p-value)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Fold change track
        fold_changes = coverage / dist_params['mean'] if dist_params['mean'] > 0 else coverage
        ax3.plot(positions, fold_changes, 'orange', alpha=0.6, linewidth=0.5)
        ax3.axhline(y=2.0, color='red', linestyle='--', alpha=0.5, 
                   label='2-fold')
        ax3.set_xlabel(f'Position on {chromosome}')
        ax3.set_ylabel('Fold Change')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def close(self):
        """Close BAM file"""
        self.bam.close()


def main():
    """Main entry point for circontrack-peaks"""
    parser = argparse.ArgumentParser(
        description="CircONTrack Coverage Peaks - Statistical peak calling for circular DNA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This module identifies statistically significant coverage peaks using negative binomial
distribution modeling, appropriate for overdispersed ONT sequencing data.

Examples:
  # Basic peak calling
  circontrack-peaks sample.bam -o peaks.bed
  
  # With reference and specific parameters
  circontrack-peaks sample.bam -r reference.fa \\
    --window-size 500 \\
    --fdr 0.01 \\
    --min-fold 3.0 \\
    -o significant_peaks.bed
  
  # Analyze specific chromosomes with plots
  circontrack-peaks sample.bam -r reference.fa \\
    --chromosomes chr1 chr2 chr3 \\
    --plot \\
    -o peaks.bed
        """
    )
    
    # Required arguments
    parser.add_argument('bam_file', 
                       help='Indexed BAM file')
    
    # Optional arguments
    parser.add_argument('-r', '--reference', 
                       help='Reference FASTA file (for chromosome lengths)')
    parser.add_argument('-o', '--output', default='coverage_peaks.bed',
                       help='Output BED file (default: coverage_peaks.bed)')
    
    # Peak calling parameters
    parser.add_argument('--window-size', type=int, default=1000,
                       help='Window size in bp (default: 1000)')
    parser.add_argument('--step-size', type=int,
                       help='Step size in bp (default: window-size/2)')
    parser.add_argument('--fdr', type=float, default=0.05,
                       help='False discovery rate threshold (default: 0.05)')
    parser.add_argument('--min-fold', type=float, default=2.0,
                       help='Minimum fold change (default: 2.0)')
    parser.add_argument('--min-coverage', type=float, default=5,
                       help='Minimum coverage in peak (default: 5)')
    parser.add_argument('--merge-distance', type=int, default=1000,
                       help='Merge peaks within distance (default: 1000)')
    
    # Analysis options
    parser.add_argument('--chromosomes', nargs='+',
                       help='Specific chromosomes to analyze')
    parser.add_argument('--plot', action='store_true',
                       help='Generate coverage plots')
    parser.add_argument('--plot-dir', default='peak_plots',
                       help='Directory for plots (default: peak_plots)')
    
    args = parser.parse_args()
    
    # Check BAM file
    if not Path(args.bam_file).exists():
        print(f"Error: BAM file not found: {args.bam_file}")
        sys.exit(1)
    
    # Check BAM index
    if not Path(f"{args.bam_file}.bai").exists():
        print("Indexing BAM file...")
        pysam.index(args.bam_file)
    
    print("CircONTrack Coverage Peaks")
    print("=" * 50)
    
    # Initialize peak caller
    peak_caller = CoveragePeakCaller(
        args.bam_file,
        args.reference,
        args.window_size,
        args.step_size
    )
    
    # Create plot directory if needed
    if args.plot:
        plot_dir = Path(args.plot_dir)
        plot_dir.mkdir(exist_ok=True)
    
    # Call peaks
    peaks = peak_caller.call_peaks_genome_wide(
        chromosomes=args.chromosomes,
        fdr_threshold=args.fdr,
        min_fold_change=args.min_fold,
        min_coverage=args.min_coverage,
        merge_distance=args.merge_distance,
        plot=args.plot
    )
    
    # Save plots if requested
    if args.plot and HAS_PLOTTING:
        for chrom, chrom_peaks in peaks.items():
            if chrom_peaks:
                # Recalculate windows for plotting
                windows = peak_caller.calculate_coverage_windows(chrom)
                coverage_values = [w['coverage'] for w in windows]
                dist_params = peak_caller.fit_negative_binomial(coverage_values)
                
                plot_file = plot_dir / f"{chrom}_coverage_peaks.png"
                peak_caller.plot_coverage_peaks(
                    chrom, windows, chrom_peaks, dist_params, plot_file
                )
    
    # Write output
    peak_caller.write_bed(peaks, args.output)
    
    # Summary statistics
    total_peaks = sum(len(p) for p in peaks.values())
    print(f"\nSummary:")
    print(f"  Total peaks: {total_peaks}")
    if peaks:
        all_peaks_flat = [p for chrom_peaks in peaks.values() for p in chrom_peaks]
        avg_fold = np.mean([p['fold_change'] for p in all_peaks_flat])
        avg_pval = np.mean([p['adjusted_pvalue'] for p in all_peaks_flat])
        print(f"  Average fold change: {avg_fold:.2f}")
        print(f"  Average adjusted p-value: {avg_pval:.2e}")
    
    # Close
    peak_caller.close()
    
    print(f"\nResults written to: {args.output}")
    if args.plot:
        print(f"Plots saved to: {args.plot_dir}/")


if __name__ == "__main__":
    main()