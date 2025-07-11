#!/usr/bin/env python3
"""
ONT-Optimized Circular DNA Detection via Multi-Modal Analysis
Main orchestrator that combines coverage patterns, junction detection, and split-read analysis
"""
import sys
import time
import logging
import argparse
from tqdm import tqdm
import pysam
import numpy as np

# Import detection modules
from .coverage_analyzer import CoverageAnalyzer
from .junction_detector import JunctionDetector
from .split_read_analyzer import SplitReadAnalyzer
from .confidence_scorer import ConfidenceScorer, MultiMethodIntegrator
from .utils import filter_candidates_by_confidence, calculate_gc_content


class CircularDNADetector:
    """Main detector that orchestrates multi-modal circular DNA detection"""
    
    def __init__(self, min_fold_enrichment=1.5, min_coverage=5, min_length=200, 
                 max_length=100000, min_confidence=0.3, verbose=True, log_level='INFO'):
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.min_confidence = min_confidence
        self.verbose = verbose
        
        # Setup logging
        self.setup_logging(log_level)
        
        # Initialize detection modules
        self.coverage_analyzer = CoverageAnalyzer(
            min_fold_enrichment=min_fold_enrichment,
            min_coverage=min_coverage,
            verbose=verbose
        )
        self.junction_detector = JunctionDetector(verbose=verbose)
        self.split_analyzer = SplitReadAnalyzer(verbose=verbose)
        
        # Initialize scoring and integration modules
        self.confidence_scorer = ConfidenceScorer()
        self.integrator = MultiMethodIntegrator()
        
    def setup_logging(self, log_level):
        """Setup logging configuration"""
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)
        
    def detect_circular_dna(self, bam_file, reference_file, output_file=None, 
                          chromosome=None):
        """
        Main detection pipeline
        """
        start_time = time.time()
        self.logger.info("Starting CircDNA detection pipeline")
        
        try:
            # Validate inputs
            self._validate_inputs(bam_file, reference_file)
            
            # Run individual detection methods
            self.logger.info("Running coverage analysis...")
            coverage_candidates = self.coverage_analyzer.detect_coverage_patterns(
                bam_file, chromosome
            )
            
            self.logger.info("Running junction detection...")
            junction_candidates = self.junction_detector.detect_junctions(
                bam_file, chromosome
            )
            
            self.logger.info("Running split-read analysis...")
            split_candidates = self.split_analyzer.analyze_split_reads(
                bam_file, chromosome
            )
            
            # Integrate results from all methods
            self.logger.info("Integrating results from all detection methods...")
            integrated_candidates = self.integrator.integrate_candidates(
                coverage_candidates, junction_candidates, split_candidates
            )
            
            # Apply confidence scoring to all candidates
            self.logger.info("Calculating confidence scores...")
            scored_candidates = self._apply_confidence_scoring(integrated_candidates)
            
            # Calculate GC content for all candidates
            self.logger.info("Calculating GC content...")
            candidates_with_gc = self._calculate_gc_content(scored_candidates, reference_file)
            
            # Filter by confidence threshold
            filtered_candidates = filter_candidates_by_confidence(
                candidates_with_gc, self.min_confidence
            )
            
            # Sort by confidence score (descending)
            final_candidates = sorted(
                filtered_candidates, 
                key=lambda x: x.confidence_score, 
                reverse=True
            )
            
            # Write output if requested
            if output_file:
                self._write_output(final_candidates, output_file)
            
            # Print summary
            self._print_summary(
                coverage_candidates, junction_candidates, split_candidates,
                integrated_candidates, final_candidates, start_time
            )
            
            return final_candidates
            
        except Exception as e:
            self.logger.error(f"Error in detection pipeline: {str(e)}")
            raise
    
    def _validate_inputs(self, bam_file, reference_file):
        """Validate input files"""
        # Check BAM file
        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                if not bam.has_index():
                    self.logger.warning("BAM file is not indexed, this will be slow")
                self.logger.info(f"✓ BAM file validated: {bam_file}")
        except Exception as e:
            raise ValueError(f"Invalid BAM file: {str(e)}")
        
        # Check reference file
        try:
            with pysam.FastaFile(reference_file) as ref:
                num_contigs = len(ref.references)
                self.logger.info(f"✓ Reference file validated: {reference_file} ({num_contigs} contigs)")
        except Exception as e:
            raise ValueError(f"Invalid reference file: {str(e)}")
    
    def _apply_confidence_scoring(self, candidates):
        """Apply confidence scoring to all candidates"""
        scored_candidates = []
        
        if self.verbose and len(candidates) > 100:
            progress = tqdm(candidates, desc="Scoring candidates", file=sys.stdout)
        else:
            progress = candidates
        
        for candidate in progress:
            # Calculate confidence score
            candidate.confidence_score = self.confidence_scorer.calculate_confidence(candidate)
            scored_candidates.append(candidate)
        
        if hasattr(progress, 'close'):
            progress.close()
        
        return scored_candidates
    
    def _calculate_gc_content(self, candidates, reference_file):
        """Calculate GC content for all candidates"""
        try:
            with pysam.FastaFile(reference_file) as ref:
                if self.verbose and len(candidates) > 100:
                    progress = tqdm(candidates, desc="Calculating GC content", file=sys.stdout)
                else:
                    progress = candidates
                
                for candidate in progress:
                    try:
                        # Fetch sequence for this candidate
                        sequence = ref.fetch(
                            candidate.chromosome, 
                            candidate.start, 
                            candidate.end
                        )
                        
                        # Calculate GC content
                        candidate.gc_content = calculate_gc_content(sequence)
                        
                    except Exception as e:
                        self.logger.warning(
                            f"Could not calculate GC content for "
                            f"{candidate.chromosome}:{candidate.start}-{candidate.end}: {str(e)}"
                        )
                        candidate.gc_content = None
                
                if hasattr(progress, 'close'):
                    progress.close()
                    
        except Exception as e:
            self.logger.error(f"Error accessing reference file for GC calculation: {str(e)}")
            # Return candidates without GC content if reference file has issues
        
        return candidates
    
    def _write_output(self, candidates, output_file):
        """Write results to BED file with confidence scores and GC content"""
        self.logger.info(f"Writing results to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("# CircDNA Detection Results\n")
            f.write("# chr\tstart\tend\tname\tconfidence\tstrand\tmethod\tlength\tgc_content\n")
            
            for i, candidate in enumerate(candidates):
                gc_str = f"{candidate.gc_content:.3f}" if candidate.gc_content is not None else "NA"
                f.write(f"{candidate.chromosome}\t{candidate.start}\t{candidate.end}\t"
                       f"circDNA_{i+1}\t{candidate.confidence_score:.3f}\t.\t"
                       f"{candidate.detection_method}\t{candidate.length}\t{gc_str}\n")
    
    def _print_summary(self, coverage_cands, junction_cands, split_cands, 
                      integrated_cands, final_cands, start_time):
        """Print comprehensive summary"""
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*60)
        print("CIRCULAR DNA DETECTION SUMMARY")
        print("="*60)
        
        print(f"Raw candidates found:")
        print(f"  Coverage-based: {len(coverage_cands)}")
        print(f"  Junction-based: {len(junction_cands)}")
        print(f"  Split-read based: {len(split_cands)}")
        print(f"  Total raw: {len(coverage_cands) + len(junction_cands) + len(split_cands)}")
        
        print(f"\nAfter multi-method integration: {len(integrated_cands)}")
        print(f"After confidence filtering (≥{self.min_confidence:.2f}): {len(final_cands)}")
        
        if final_cands:
            print(f"\nTop 5 candidates by confidence:")
            for i, candidate in enumerate(final_cands[:5]):
                print(f"  {i+1}. {candidate.chromosome}:{candidate.start}-{candidate.end}")
                print(f"     Confidence: {candidate.confidence_score:.3f}")
                print(f"     Method: {candidate.detection_method}")
                print(f"     Length: {candidate.length:,} bp")
                if candidate.gc_content is not None:
                    print(f"     GC content: {candidate.gc_content:.1%}")
                
                # Show contributing evidence
                evidence = []
                if candidate.fold_enrichment:
                    evidence.append(f"fold={candidate.fold_enrichment:.1f}")
                if candidate.junction_support:
                    evidence.append(f"junctions={candidate.junction_support}")
                if candidate.split_support:
                    evidence.append(f"splits={candidate.split_support}")
                
                if evidence:
                    print(f"     Evidence: {', '.join(evidence)}")
                print()
        
        print(f"\nAnalysis completed in: {elapsed_time:.2f} seconds")
        print("="*60)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Detect circular DNA from ONT sequencing data using multi-modal analysis'
    )
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('reference_file', help='Reference FASTA file')
    parser.add_argument('-o', '--output', help='Output BED file', 
                       default='circular_dna_results.bed')
    parser.add_argument('-c', '--chromosome', help='Analyze specific chromosome only')
    parser.add_argument('-q', '--quiet', action='store_true', 
                       help='Disable verbose output (quiet mode)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                       default='INFO', help='Set logging level')
    
    # Detection parameters
    parser.add_argument('--min-fold-enrichment', type=float, default=1.5,
                       help='Minimum fold enrichment for coverage-based detection')
    parser.add_argument('--min-coverage', type=int, default=5,
                       help='Minimum coverage threshold')
    parser.add_argument('--min-length', type=int, default=200,
                       help='Minimum circular DNA length')
    parser.add_argument('--max-length', type=int, default=100000,
                       help='Maximum circular DNA length')
    parser.add_argument('--min-confidence', type=float, default=0.3,
                       help='Minimum confidence score threshold')

    args = parser.parse_args()
    
    # Create detector
    detector = CircularDNADetector(
        min_fold_enrichment=args.min_fold_enrichment,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        max_length=args.max_length,
        min_confidence=args.min_confidence,
        verbose=not args.quiet,
        log_level=args.log_level
    )
    
    # Run detection
    results = detector.detect_circular_dna(
        args.bam_file,
        args.reference_file,
        args.output,
        args.chromosome
    )
    
    print(f"\nAnalysis complete! Results written to {args.output}")
    print(f"Found {len(results)} high-confidence circular DNA candidates")


if __name__ == "__main__":
    main()