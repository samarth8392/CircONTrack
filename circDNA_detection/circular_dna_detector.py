#!/usr/bin/env python3
"""
ONT-Optimized Circular DNA Detection via Multi-Modal Analysis
Combines coverage patterns, junction detection, and split-read analysis
"""
import sys
import time
import logging
from tqdm import tqdm
import pysam
import numpy as np

class CircularDNADetector:
    """Enhanced detector with comprehensive progress tracking and verbose output"""
    
    def __init__(self, min_fold_enrichment=1.5, min_coverage=5, min_length=200, 
                 max_length=100000, verbose=True, log_level='INFO'):  # Changed default to True
        self.min_fold_enrichment = min_fold_enrichment
        self.min_coverage = min_coverage
        self.min_length = min_length
        self.max_length = max_length
        self.verbose = verbose
        
        # Setup logging
        self.setup_logging(log_level)
        
        # Progress tracking
        self.total_steps = 0
        self.current_step = 0
        self.step_descriptions = []
        
    def setup_logging(self, log_level):
        """Setup logging configuration"""
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)
        
    def log_step(self, message, level='INFO'):
        """Log a step with progress information"""
        if self.verbose:
            progress = f"[{self.current_step}/{self.total_steps}]" if self.total_steps > 0 else ""
            getattr(self.logger, level.lower())(f"{progress} {message}")
            
    def detect_circular_dna(self, bam_file, reference_file, output_file=None, 
                          chromosome=None):
        """
        Main detection pipeline with comprehensive progress tracking
        """
        start_time = time.time()
        self.log_step("Starting CircDNA detection pipeline", "INFO")
        
        # Initialize progress tracking
        self._initialize_progress_tracking(bam_file, chromosome)
        
        try:
            # Phase 1: File validation and setup
            self._validate_inputs(bam_file, reference_file)
            
            # Phase 2: BAM file analysis
            bam_stats = self._analyze_bam_file(bam_file)
            
            # Phase 3: Process each chromosome
            candidates = self._process_chromosomes(bam_file, reference_file, chromosome)
            
            # Phase 4: Final processing and output
            final_results = self._finalize_results(candidates, output_file)
            
            # Summary
            self._print_summary(final_results, start_time)
            
            return final_results
            
        except Exception as e:
            self.logger.error(f"Error in detection pipeline: {str(e)}")
            raise
    
    def _initialize_progress_tracking(self, bam_file, chromosome):
        """Initialize progress tracking based on input"""
        self.log_step("Initializing progress tracking...")
        
        # Estimate total steps based on input
        if chromosome:
            self.total_steps = 8  # Single chromosome
        else:
            # Estimate based on BAM file
            try:
                with pysam.AlignmentFile(bam_file, 'rb') as bam:
                    num_chromosomes = len(bam.references)
                    self.total_steps = 4 + (num_chromosomes * 4)  # Setup + per-chr steps
            except:
                self.total_steps = 50  # Default estimate
                
        self.current_step = 0
        self.log_step(f"Estimated {self.total_steps} total steps")
    
    def _validate_inputs(self, bam_file, reference_file):
        """Validate input files with progress updates"""
        self.current_step += 1
        self.log_step(f"Validating input files...")
        
        # Check BAM file
        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                if not bam.has_index():
                    self.log_step("WARNING: BAM file is not indexed, this will be slow", "WARNING")
                self.log_step(f"✓ BAM file validated: {bam_file}")
        except Exception as e:
            raise ValueError(f"Invalid BAM file: {str(e)}")
        
        # Check reference file
        try:
            with pysam.FastaFile(reference_file) as ref:
                num_contigs = len(ref.references)
                self.log_step(f"✓ Reference file validated: {reference_file} ({num_contigs} contigs)")
        except Exception as e:
            raise ValueError(f"Invalid reference file: {str(e)}")
    
    def _analyze_bam_file(self, bam_file):
        """Analyze BAM file statistics with progress"""
        self.current_step += 1
        self.log_step("Analyzing BAM file statistics...")
        
        stats = {
            'total_reads': 0,
            'mapped_reads': 0,
            'chromosomes': [],
            'avg_coverage': 0
        }
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            # Get chromosome information
            stats['chromosomes'] = list(bam.references)
            self.log_step(f"Found {len(stats['chromosomes'])} chromosomes/contigs")
            
            if self.verbose:
                # Quick read count (sample-based for speed)
                self.log_step("Sampling reads for statistics...")
                sample_size = 10000
                read_count = 0
                mapped_count = 0
                
                # Progress bar for read sampling
                sample_progress = tqdm(desc="Sampling reads", 
                                     total=sample_size, 
                                     disable=not self.verbose, 
                                     file=sys.stdout)
                
                for i, read in enumerate(bam.fetch()):
                    if i >= sample_size:
                        break
                    read_count += 1
                    if not read.is_unmapped:
                        mapped_count += 1
                    sample_progress.update(1)
                
                sample_progress.close()
                
                if read_count > 0:
                    stats['mapped_percentage'] = (mapped_count / read_count) * 100
                    self.log_step(f"Sample stats: {mapped_count}/{read_count} reads mapped "
                                f"({stats['mapped_percentage']:.1f}%)")
        
        return stats
    
    def _process_chromosomes(self, bam_file, reference_file, target_chromosome=None):
        """Process chromosomes with detailed progress tracking"""
        self.current_step += 1
        self.log_step("Starting chromosome-by-chromosome analysis...")
        
        all_candidates = []
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            chromosomes = [target_chromosome] if target_chromosome else bam.references
            
            # Create chromosome progress bar - always show when verbose is true
            chr_progress = tqdm(chromosomes, 
                              desc="Processing chromosomes", 
                              disable=not self.verbose, 
                              file=sys.stdout,
                              position=0,
                              leave=True)
            
            for chr_name in chr_progress:
                chr_progress.set_description(f"Processing {chr_name}")
                
                try:
                    # Get chromosome length
                    chr_length = bam.get_reference_length(chr_name)
                    self.log_step(f"Analyzing chromosome {chr_name} (length: {chr_length:,} bp)")
                    
                    # Phase 1: Coverage analysis
                    coverage_candidates = self._analyze_coverage(bam, chr_name, chr_length)
                    
                    # Phase 2: Junction detection
                    junction_candidates = self._detect_junctions(bam, chr_name)
                    
                    # Phase 3: Split-read analysis
                    split_candidates = self._analyze_split_reads(bam, chr_name)
                    
                    # Phase 4: Integrate results for this chromosome
                    chr_candidates = self._integrate_chromosome_results(
                        chr_name, coverage_candidates, junction_candidates, split_candidates
                    )
                    
                    all_candidates.extend(chr_candidates)
                    
                    self.log_step(f"Chromosome {chr_name} complete: {len(chr_candidates)} candidates found")
                    
                except Exception as e:
                    self.log_step(f"Error processing chromosome {chr_name}: {str(e)}", "ERROR")
                    continue
            
            chr_progress.close()
        
        return all_candidates
    
    def _analyze_coverage(self, bam, chromosome, chr_length):
        """Analyze coverage patterns with progress tracking"""
        self.current_step += 1
        self.log_step(f"  → Phase 1: Coverage analysis for {chromosome}")
        
        candidates = []
        window_size = 1000  # 1kb windows
        num_windows = chr_length // window_size + 1
        
        # Progress bar for coverage analysis - always show when verbose
        coverage_progress = tqdm(range(0, chr_length, window_size), 
                               desc=f"Coverage analysis ({chromosome})", 
                               disable=not self.verbose,
                               file=sys.stdout,
                               position=1,
                               leave=False)
        
        for start in coverage_progress:
            end = min(start + window_size, chr_length)
            
            try:
                # Calculate coverage for this window
                coverage = bam.count_coverage(chromosome, start, end)
                total_coverage = sum(sum(base_cov) for base_cov in coverage)
                avg_coverage = total_coverage / (end - start) if end > start else 0
                
                # Check if this window meets criteria
                if avg_coverage >= self.min_coverage:
                    # Additional fold-enrichment check would go here
                    candidates.append({
                        'chr': chromosome,
                        'start': start,
                        'end': end,
                        'coverage': avg_coverage,
                        'method': 'coverage'
                    })
                    
            except Exception as e:
                if self.verbose:
                    self.log_step(f"    Error in coverage window {start}-{end}: {str(e)}", "DEBUG")
                continue
        
        coverage_progress.close()
        self.log_step(f"  → Coverage analysis complete: {len(candidates)} regions with elevated coverage")
        return candidates
    
    def _detect_junctions(self, bam, chromosome):
        """Detect junction signatures with progress tracking"""
        self.current_step += 1
        self.log_step(f"  → Phase 2: Junction detection for {chromosome}")
        
        candidates = []
        junction_count = 0
        
        # Get total reads for this chromosome for progress tracking
        total_reads = bam.count(chromosome)
        
        # Progress tracking for junction detection
        junction_progress = tqdm(bam.fetch(chromosome), 
                               desc=f"Junction detection ({chromosome})", 
                               total=total_reads,
                               disable=not self.verbose,
                               file=sys.stdout,
                               position=1,
                               leave=False)
        
        try:
            for read in junction_progress:
                # Junction detection logic
                if self._is_junction_read(read):
                    junction_count += 1
                    candidates.append({
                        'chr': chromosome,
                        'start': read.reference_start,
                        'end': read.reference_end,
                        'method': 'junction',
                        'read_name': read.query_name
                    })
                    
                    # Update progress description with current count
                    if junction_count % 100 == 0:
                        junction_progress.set_description(f"Junction detection ({chromosome}) - Found: {junction_count}")
                    
        except Exception as e:
            self.log_step(f"Error in junction detection: {str(e)}", "ERROR")
        finally:
            junction_progress.close()
        
        self.log_step(f"  → Junction detection complete: {len(candidates)} junction signatures found")
        return candidates
    
    def _analyze_split_reads(self, bam, chromosome):
        """Analyze split reads with progress tracking"""
        self.current_step += 1
        self.log_step(f"  → Phase 3: Split-read analysis for {chromosome}")
        
        candidates = []
        split_count = 0
        
        # Get total reads for progress tracking
        total_reads = bam.count(chromosome)
        
        # Progress bar for split-read analysis
        split_progress = tqdm(bam.fetch(chromosome), 
                            desc=f"Split-read analysis ({chromosome})", 
                            total=total_reads,
                            disable=not self.verbose,
                            file=sys.stdout,
                            position=1,
                            leave=False)
        
        try:
            for read in split_progress:
                if self._is_split_read(read):
                    split_count += 1
                    candidates.append({
                        'chr': chromosome,
                        'start': read.reference_start,
                        'end': read.reference_end,
                        'method': 'split_read',
                        'read_name': read.query_name
                    })
                    
                    # Update progress description with current count
                    if split_count % 100 == 0:
                        split_progress.set_description(f"Split-read analysis ({chromosome}) - Found: {split_count}")
                    
        except Exception as e:
            self.log_step(f"Error in split-read analysis: {str(e)}", "ERROR")
        finally:
            split_progress.close()
        
        self.log_step(f"  → Split-read analysis complete: {len(candidates)} split-read signatures found")
        return candidates
    
    def _integrate_chromosome_results(self, chromosome, coverage_cands, junction_cands, split_cands):
        """Integrate results from all three methods"""
        self.current_step += 1
        self.log_step(f"  → Phase 4: Integrating results for {chromosome}")
        
        # Simple integration logic (you'd make this more sophisticated)
        all_candidates = coverage_cands + junction_cands + split_cands
        
        # Filter by length requirements with progress bar
        filtered_candidates = []
        
        if self.verbose and len(all_candidates) > 1000:
            # Show progress bar for filtering if there are many candidates
            filter_progress = tqdm(all_candidates, 
                                 desc=f"Filtering candidates ({chromosome})", 
                                 disable=not self.verbose,
                                 file=sys.stdout,
                                 position=1,
                                 leave=False)
        else:
            filter_progress = all_candidates
        
        for cand in filter_progress:
            length = cand['end'] - cand['start']
            if self.min_length <= length <= self.max_length:
                cand['length'] = length
                filtered_candidates.append(cand)
        
        if hasattr(filter_progress, 'close'):
            filter_progress.close()
        
        self.log_step(f"  → Integration complete: {len(filtered_candidates)} candidates after filtering")
        return filtered_candidates
    
    def _finalize_results(self, candidates, output_file):
        """Finalize results and write output"""
        self.current_step += 1
        self.log_step("Finalizing results...")
        
        # Sort candidates by confidence/evidence with progress bar
        if self.verbose and len(candidates) > 1000:
            self.log_step("Sorting candidates by length...")
            # For very large datasets, we might want to show sorting progress
            
        sorted_candidates = sorted(candidates, key=lambda x: x.get('length', 0), reverse=True)
        
        if output_file:
            self._write_output(sorted_candidates, output_file)
        
        self.log_step(f"Results finalized: {len(sorted_candidates)} total candidates")
        return sorted_candidates
    
    def _write_output(self, candidates, output_file):
        """Write results to BED file"""
        self.log_step(f"Writing results to {output_file}")
        
        # Progress bar for writing output if many candidates
        if self.verbose and len(candidates) > 1000:
            write_progress = tqdm(candidates, 
                                desc="Writing output", 
                                disable=not self.verbose,
                                file=sys.stdout)
        else:
            write_progress = candidates
        
        with open(output_file, 'w') as f:
            f.write("# CircDNA Detection Results\n")
            f.write("# chr\tstart\tend\tname\tscore\tstrand\tmethod\tlength\n")
            
            for i, cand in enumerate(write_progress):
                f.write(f"{cand['chr']}\t{cand['start']}\t{cand['end']}\t"
                       f"circDNA_{i+1}\t{cand.get('coverage', 0):.1f}\t.\t"
                       f"{cand['method']}\t{cand.get('length', 0)}\n")
        
        if hasattr(write_progress, 'close'):
            write_progress.close()
    
    def _print_summary(self, results, start_time):
        """Print comprehensive summary"""
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*60)
        print("CIRCULAR DNA DETECTION SUMMARY")
        print("="*60)
        print(f"Total candidates found: {len(results)}")
        print(f"Analysis completed in: {elapsed_time:.2f} seconds")
        
        if results:
            # Method breakdown
            method_counts = {}
            for result in results:
                method = result['method']
                method_counts[method] = method_counts.get(method, 0) + 1
            
            print("\nCandidates by detection method:")
            for method, count in method_counts.items():
                print(f"  {method}: {count}")
            
            # Size distribution
            lengths = [r.get('length', 0) for r in results]
            if lengths:
                print(f"\nSize distribution:")
                print(f"  Min length: {min(lengths):,} bp")
                print(f"  Max length: {max(lengths):,} bp")
                print(f"  Mean length: {np.mean(lengths):,.0f} bp")
        
        print("="*60)
    
    def _is_junction_read(self, read):
        """Check if read shows junction signature"""
        # Simplified junction detection logic
        return (not read.is_unmapped and 
                read.mapping_quality >= 20 and
                read.cigarstring and 'S' in read.cigarstring)
    
    def _is_split_read(self, read):
        """Check if read is split alignment"""
        # Simplified split-read detection
        return (read.has_tag('SA') or  # Supplementary alignment
                (read.cigarstring and read.cigarstring.count('S') >= 2))


def main():
    """Enhanced main function with argument parsing"""
    import argparse
    
    parser = argparse.ArgumentParser(description='CircDNA Detection with Enhanced Progress Tracking')
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('reference_file', help='Reference FASTA file')
    parser.add_argument('-o', '--output', help='Output BED file', default='circular_dna_verbose.bed')
    parser.add_argument('-c', '--chromosome', help='Analyze specific chromosome only')
    parser.add_argument('-q', '--quiet', action='store_true', help='Disable verbose output (quiet mode)')  # Changed to quiet mode
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                       default='INFO', help='Set logging level')
    parser.add_argument('--min-fold-enrichment', type=float, default=1.5)
    parser.add_argument('--min-coverage', type=int, default=5)
    parser.add_argument('--min-length', type=int, default=200)
    parser.add_argument('--max-length', type=int, default=100000)

    args = parser.parse_args()
    
    # Create detector with enhanced options - verbose is now default unless quiet is specified
    detector = CircularDNADetector(
        min_fold_enrichment=args.min_fold_enrichment,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        max_length=args.max_length,
        verbose=not args.quiet,  # Verbose is default, disabled only with --quiet
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


if __name__ == "__main__":
    main()