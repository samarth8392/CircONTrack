#!/usr/bin/env python3
"""
Comprehensive Peak Filtering Pipeline for CircONTrack
Removes repetitive element artifacts and validates true eccDNA peaks
"""

import pysam
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

class PeakFilterPipeline:
    """Filter CircONTrack peaks to remove artifacts"""
    
    def __init__(self, peak_file, bam_file, output_prefix="filtered"):
        self.peak_file = peak_file
        self.bam_file = bam_file
        self.output_prefix = output_prefix
        
        # Load peaks
        self.peaks = pd.read_csv(peak_file, sep='\t',
                                 names=['chr', 'start', 'end', 'name', 'score',
                                        'strand', 'coverage', 'fold_change',
                                        'pvalue', 'adjusted_pvalue', 'read_count'],
                                 comment='#')
        
        self.bam = pysam.AlignmentFile(bam_file, 'rb')
        self.filtered_peaks = None
        self.artifact_peaks = None
        
    def analyze_peak_reads(self, chrom, start, end, max_reads=1000):
        """
        Analyze read characteristics in a peak region
        
        Returns dict with quality metrics
        """
        metrics = {
            'total_reads': 0,
            'soft_clipped': 0,
            'high_soft_clip': 0,  # >30% of read soft-clipped
            'extreme_soft_clip': 0,  # >50% of read soft-clipped
            'low_mapq': 0,  # MAPQ < 20
            'very_low_mapq': 0,  # MAPQ < 10
            'mean_mapq': 0,
            'unique_starts': set(),
            'unique_ends': set(),
            'complex_cigars': 0,  # More than 3 CIGAR operations
            'max_soft_clip_pct': 0,
            'median_soft_clip_pct': []
        }
        
        mapq_scores = []
        
        for i, read in enumerate(self.bam.fetch(chrom, start, end)):
            if i >= max_reads:
                break
                
            metrics['total_reads'] += 1
            metrics['unique_starts'].add(read.reference_start)
            metrics['unique_ends'].add(read.reference_end)
            
            # MAPQ analysis
            mapq_scores.append(read.mapping_quality)
            if read.mapping_quality < 20:
                metrics['low_mapq'] += 1
            if read.mapping_quality < 10:
                metrics['very_low_mapq'] += 1
            
            # CIGAR analysis
            if read.cigartuples:
                # Calculate soft-clip percentage
                soft_clip = sum(length for op, length in read.cigartuples if op == 4)
                total_length = sum(length for op, length in read.cigartuples if op in [0, 1, 4])
                
                if total_length > 0:
                    soft_clip_pct = (soft_clip / total_length) * 100
                    metrics['median_soft_clip_pct'].append(soft_clip_pct)
                    
                    if soft_clip_pct > 30:
                        metrics['high_soft_clip'] += 1
                    if soft_clip_pct > 50:
                        metrics['extreme_soft_clip'] += 1
                    if soft_clip_pct > metrics['max_soft_clip_pct']:
                        metrics['max_soft_clip_pct'] = soft_clip_pct
                
                # Check CIGAR complexity
                if len(read.cigartuples) > 3:
                    metrics['complex_cigars'] += 1
                    
                # Check for specific bad patterns
                cigar_str = read.cigarstring if read.cigarstring else ""
                if 'S' in cigar_str:
                    metrics['soft_clipped'] += 1
        
        # Calculate summary statistics
        if mapq_scores:
            metrics['mean_mapq'] = np.mean(mapq_scores)
        
        if metrics['median_soft_clip_pct']:
            metrics['median_soft_clip_pct'] = np.median(metrics['median_soft_clip_pct'])
        else:
            metrics['median_soft_clip_pct'] = 0
            
        metrics['unique_starts'] = len(metrics['unique_starts'])
        metrics['unique_ends'] = len(metrics['unique_ends'])
        
        # Calculate quality flags
        metrics['soft_clip_rate'] = metrics['soft_clipped'] / metrics['total_reads'] if metrics['total_reads'] > 0 else 0
        metrics['high_soft_clip_rate'] = metrics['high_soft_clip'] / metrics['total_reads'] if metrics['total_reads'] > 0 else 0
        metrics['extreme_soft_clip_rate'] = metrics['extreme_soft_clip'] / metrics['total_reads'] if metrics['total_reads'] > 0 else 0
        metrics['low_mapq_rate'] = metrics['low_mapq'] / metrics['total_reads'] if metrics['total_reads'] > 0 else 0
        metrics['position_diversity'] = min(metrics['unique_starts'], metrics['unique_ends']) / metrics['total_reads'] if metrics['total_reads'] > 0 else 0
        
        return metrics
    
    def classify_peaks(self):
        """
        Classify peaks as likely eccDNA vs artifacts
        """
        print("Analyzing all peaks for artifacts...")
        print("-" * 60)
        
        classifications = []
        
        for idx, peak in self.peaks.iterrows():
            if idx % 10 == 0:
                print(f"  Processing peak {idx+1}/{len(self.peaks)}...")
            
            # Analyze reads in this peak
            metrics = self.analyze_peak_reads(peak['chr'], peak['start'], peak['end'])
            
            # Classification rules
            is_artifact = False
            artifact_reasons = []
            confidence_score = 100  # Start with 100, deduct for issues
            
            # Rule 1: Extreme fold change
            if peak['fold_change'] > 100:
                is_artifact = True
                artifact_reasons.append(f"Extreme fold change ({peak['fold_change']:.0f}x)")
                confidence_score -= 50
            elif peak['fold_change'] > 50:
                artifact_reasons.append(f"High fold change ({peak['fold_change']:.0f}x)")
                confidence_score -= 25
            
            # Rule 2: Soft clipping
            if metrics['extreme_soft_clip_rate'] > 0.3:
                is_artifact = True
                artifact_reasons.append(f"Extreme soft-clipping ({metrics['extreme_soft_clip_rate']:.1%})")
                confidence_score -= 40
            elif metrics['high_soft_clip_rate'] > 0.3:
                artifact_reasons.append(f"High soft-clipping ({metrics['high_soft_clip_rate']:.1%})")
                confidence_score -= 20
            
            # Rule 3: Mapping quality
            if metrics['mean_mapq'] < 10:
                is_artifact = True
                artifact_reasons.append(f"Very low MAPQ ({metrics['mean_mapq']:.1f})")
                confidence_score -= 50
            elif metrics['mean_mapq'] < 20:
                artifact_reasons.append(f"Low MAPQ ({metrics['mean_mapq']:.1f})")
                confidence_score -= 25
            
            # Rule 4: Position diversity
            if metrics['position_diversity'] < 0.1:
                is_artifact = True
                artifact_reasons.append("Low position diversity (stacked reads)")
                confidence_score -= 30
            
            # Rule 5: Read count (coverage)
            if metrics['total_reads'] > 10000:
                is_artifact = True
                artifact_reasons.append(f"Extreme coverage ({metrics['total_reads']} reads)")
                confidence_score -= 40
            
            # Store classification
            classification = {
                'peak_name': peak['name'],
                'chr': peak['chr'],
                'start': peak['start'],
                'end': peak['end'],
                'fold_change': peak['fold_change'],
                'coverage': peak['coverage'],
                'pvalue': peak['adjusted_pvalue'],
                'is_artifact': is_artifact,
                'confidence_score': max(0, confidence_score),
                'artifact_reasons': '; '.join(artifact_reasons) if artifact_reasons else 'None',
                'mean_mapq': metrics['mean_mapq'],
                'soft_clip_rate': metrics['high_soft_clip_rate'],
                'position_diversity': metrics['position_diversity'],
                'read_count': metrics['total_reads']
            }
            
            classifications.append(classification)
        
        self.classifications = pd.DataFrame(classifications)
        
        # Separate into filtered and artifact peaks
        self.filtered_peaks = self.classifications[~self.classifications['is_artifact']]
        self.artifact_peaks = self.classifications[self.classifications['is_artifact']]
        
        print(f"\nClassification complete:")
        print(f"  Total peaks: {len(self.peaks)}")
        print(f"  Likely eccDNA: {len(self.filtered_peaks)}")
        print(f"  Artifacts: {len(self.artifact_peaks)}")
        
        return self.classifications
    
    def generate_report(self):
        """Generate filtering report with visualizations"""
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Plot 1: Fold change distribution
        ax1 = axes[0, 0]
        ax1.hist(self.filtered_peaks['fold_change'], bins=30, alpha=0.7, 
                label='Kept', color='green', edgecolor='black')
        ax1.hist(self.artifact_peaks['fold_change'], bins=30, alpha=0.7,
                label='Filtered', color='red', edgecolor='black')
        ax1.set_xlabel('Fold Change')
        ax1.set_ylabel('Count')
        ax1.set_title('Fold Change Distribution')
        ax1.legend()
        ax1.set_yscale('log')
        
        # Plot 2: MAPQ distribution
        ax2 = axes[0, 1]
        ax2.scatter(self.filtered_peaks['mean_mapq'], self.filtered_peaks['fold_change'],
                   alpha=0.6, color='green', label='Kept', s=20)
        ax2.scatter(self.artifact_peaks['mean_mapq'], self.artifact_peaks['fold_change'],
                   alpha=0.6, color='red', label='Filtered', s=20)
        ax2.set_xlabel('Mean MAPQ')
        ax2.set_ylabel('Fold Change')
        ax2.set_title('MAPQ vs Fold Change')
        ax2.axvline(x=20, color='black', linestyle='--', alpha=0.3)
        ax2.legend()
        
        # Plot 3: Soft-clipping rate
        ax3 = axes[0, 2]
        ax3.scatter(self.filtered_peaks['soft_clip_rate']*100, self.filtered_peaks['fold_change'],
                   alpha=0.6, color='green', label='Kept', s=20)
        ax3.scatter(self.artifact_peaks['soft_clip_rate']*100, self.artifact_peaks['fold_change'],
                   alpha=0.6, color='red', label='Filtered', s=20)
        ax3.set_xlabel('High Soft-clip Rate (%)')
        ax3.set_ylabel('Fold Change')
        ax3.set_title('Soft-clipping vs Fold Change')
        ax3.axvline(x=30, color='black', linestyle='--', alpha=0.3)
        ax3.legend()
        
        # Plot 4: Confidence scores
        ax4 = axes[1, 0]
        all_scores = pd.concat([self.filtered_peaks, self.artifact_peaks])
        ax4.hist(self.filtered_peaks['confidence_score'], bins=20, alpha=0.7,
                label='Kept', color='green', edgecolor='black')
        ax4.hist(self.artifact_peaks['confidence_score'], bins=20, alpha=0.7,
                label='Filtered', color='red', edgecolor='black')
        ax4.set_xlabel('Confidence Score')
        ax4.set_ylabel('Count')
        ax4.set_title('Peak Confidence Distribution')
        ax4.legend()
        
        # Plot 5: Position diversity
        ax5 = axes[1, 1]
        ax5.scatter(self.filtered_peaks['position_diversity'], self.filtered_peaks['fold_change'],
                   alpha=0.6, color='green', label='Kept', s=20)
        ax5.scatter(self.artifact_peaks['position_diversity'], self.artifact_peaks['fold_change'],
                   alpha=0.6, color='red', label='Filtered', s=20)
        ax5.set_xlabel('Position Diversity')
        ax5.set_ylabel('Fold Change')
        ax5.set_title('Read Position Diversity')
        ax5.axvline(x=0.1, color='black', linestyle='--', alpha=0.3)
        ax5.legend()
        
        # Plot 6: Summary statistics
        ax6 = axes[1, 2]
        ax6.axis('off')
        
        summary_text = f"""
        Filtering Summary
        {'='*30}
        
        Input peaks: {len(self.peaks)}
        Kept (likely eccDNA): {len(self.filtered_peaks)}
        Filtered (artifacts): {len(self.artifact_peaks)}
        
        Filter rate: {len(self.artifact_peaks)/len(self.peaks)*100:.1f}%
        
        Top artifact reasons:
        """
        
        # Count artifact reasons
        if len(self.artifact_peaks) > 0:
            reasons = []
            for r in self.artifact_peaks['artifact_reasons']:
                reasons.extend(r.split('; '))
            
            from collections import Counter
            reason_counts = Counter(reasons)
            for reason, count in reason_counts.most_common(5):
                if reason != 'None':
                    summary_text += f"\n  • {reason}: {count}"
        
        ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        plt.suptitle('CircONTrack Peak Filtering Report', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # Save figure
        plt.savefig(f"{self.output_prefix}_filtering_report.png", dpi=150, bbox_inches='tight')
        print(f"\nReport saved to {self.output_prefix}_filtering_report.png")
        
        plt.show()
    
    def save_results(self):
        """Save filtered peaks and classification results"""
        
        # Save filtered peaks in BED format
        if len(self.filtered_peaks) > 0:
            filtered_bed = self.peaks[self.peaks['name'].isin(self.filtered_peaks['peak_name'])]
            filtered_bed.to_csv(f"{self.output_prefix}_peaks.bed", 
                              sep='\t', index=False, header=False)
            print(f"Filtered peaks saved to {self.output_prefix}_peaks.bed")
        
        # Save artifact peaks
        if len(self.artifact_peaks) > 0:
            artifact_bed = self.peaks[self.peaks['name'].isin(self.artifact_peaks['peak_name'])]
            artifact_bed.to_csv(f"{self.output_prefix}_artifacts.bed",
                              sep='\t', index=False, header=False)
            print(f"Artifact peaks saved to {self.output_prefix}_artifacts.bed")
        
        # Save full classification table
        self.classifications.to_csv(f"{self.output_prefix}_classifications.tsv",
                                   sep='\t', index=False)
        print(f"Full classifications saved to {self.output_prefix}_classifications.tsv")
        
        # Save high-confidence eccDNA candidates
        high_conf = self.filtered_peaks[self.filtered_peaks['confidence_score'] >= 75]
        if len(high_conf) > 0:
            high_conf_bed = self.peaks[self.peaks['name'].isin(high_conf['peak_name'])]
            high_conf_bed.to_csv(f"{self.output_prefix}_high_confidence.bed",
                               sep='\t', index=False, header=False)
            print(f"High-confidence peaks saved to {self.output_prefix}_high_confidence.bed")
        
        # Generate summary report
        with open(f"{self.output_prefix}_summary.txt", 'w') as f:
            f.write("CircONTrack Peak Filtering Summary\n")
            f.write("="*50 + "\n\n")
            
            f.write(f"Total input peaks: {len(self.peaks)}\n")
            f.write(f"Peaks kept (likely eccDNA): {len(self.filtered_peaks)}\n")
            f.write(f"Peaks filtered (artifacts): {len(self.artifact_peaks)}\n")
            f.write(f"Filter rate: {len(self.artifact_peaks)/len(self.peaks)*100:.1f}%\n\n")
            
            if len(self.artifact_peaks) > 0:
                f.write("Top 10 filtered peaks (worst artifacts):\n")
                f.write("-"*50 + "\n")
                for _, peak in self.artifact_peaks.head(10).iterrows():
                    f.write(f"{peak['chr']}:{peak['start']}-{peak['end']}\n")
                    f.write(f"  Fold change: {peak['fold_change']:.1f}x\n")
                    f.write(f"  Reasons: {peak['artifact_reasons']}\n\n")
            
            if len(high_conf) > 0:
                f.write(f"\nHigh-confidence eccDNA candidates: {len(high_conf)}\n")
                f.write("-"*50 + "\n")
                for _, peak in high_conf.head(10).iterrows():
                    f.write(f"{peak['chr']}:{peak['start']}-{peak['end']}\n")
                    f.write(f"  Fold change: {peak['fold_change']:.1f}x\n")
                    f.write(f"  Confidence: {peak['confidence_score']:.0f}%\n\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter CircONTrack peaks to remove artifacts"
    )
    parser.add_argument('peak_file', help='Peak BED file from CircONTrack')
    parser.add_argument('bam_file', help='BAM file used for peak calling')
    parser.add_argument('-o', '--output-prefix', default='filtered',
                       help='Output file prefix')
    parser.add_argument('--max-fold', type=float, default=100,
                       help='Maximum fold change allowed')
    parser.add_argument('--min-mapq', type=float, default=20,
                       help='Minimum mean MAPQ required')
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = PeakFilterPipeline(args.peak_file, args.bam_file, args.output_prefix)
    
    # Run classification
    pipeline.classify_peaks()
    
    # Generate report
    pipeline.generate_report()
    
    # Save results
    pipeline.save_results()
    
    print("\nFiltering complete!")


if __name__ == "__main__":
    main()