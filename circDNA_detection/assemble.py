#!/usr/bin/env python3
"""
CircONTrack Assemble Module - Extract and prepare reads for circular DNA assembly
Part of the CircONTrack pipeline for circular DNA analysis
"""

import os
import sys
import pysam
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

# Try to import pandas for stats, but make it optional
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# Try to import tqdm for progress bars, but make it optional
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    # Simple progress printer if tqdm not available
    def tqdm(iterable, total=None, desc=None):
        return iterable

class CircularDNAAssembler:
    """
    Extract reads from circular DNA regions and prepare for assembly
    Part of the CircONTrack suite for circular DNA analysis
    """
    
    def __init__(self, bam_file, output_dir, threads=4, padding=1000):
        """
        Initialize the assembler module
        
        Parameters:
        -----------
        bam_file : str
            Path to input BAM file (should be same used for CircONTrack detection)
        output_dir : str
            Directory for output files
        threads : int
            Number of threads for parallel processing
        padding : int
            Extra bases to include around regions for better assembly
        """
        self.bam_file = bam_file
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.padding = padding
        
        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.reads_dir = self.output_dir / "reads"
        self.reads_dir.mkdir(exist_ok=True)
        self.info_dir = self.output_dir / "info"
        self.info_dir.mkdir(exist_ok=True)
        
        # Check if BAM index exists
        if not os.path.exists(f"{bam_file}.bai"):
            print(f"Indexing BAM file...")
            pysam.index(bam_file)
        
        print(f"CircONTrack Assemble initialized")
        print(f"Output directory: {self.output_dir}")
    
    def extract_reads_from_region(self, chrom, start, end, region_id, region_type=None):
        """
        Extract reads from a specific circular DNA region
        
        Parameters:
        -----------
        chrom : str
            Chromosome/contig name
        start : int
            Region start position
        end : int
            Region end position
        region_id : str
            Unique identifier for the region
        region_type : str
            Type from classification (e.g., 'host_eccDNA', 'viral_episome')
        
        Returns:
        --------
        dict with region info and statistics
        """
        # Add padding for better assembly context
        padded_start = max(0, start - self.padding)
        padded_end = end + self.padding
        
        # Output files
        fastq_file = self.reads_dir / f"{region_id}.fastq"
        info_file = self.info_dir / f"{region_id}.info"
        
        stats = {
            'region_id': region_id,
            'chrom': chrom,
            'start': start,
            'end': end,
            'padded_start': padded_start,
            'padded_end': padded_end,
            'length': end - start,
            'region_type': region_type or 'unknown',
            'reads_extracted': 0,
            'primary_reads': 0,
            'supplementary_reads': 0,
            'reads_spanning_junction': 0,
            'chimeric_reads': 0,
            'avg_read_length': 0,
            'coverage_depth': 0,
            'read_n50': 0
        }
        
        try:
            with pysam.AlignmentFile(self.bam_file, "rb") as bam:
                with open(fastq_file, 'w') as fq:
                    read_lengths = []
                    read_names = set()
                    
                    # Extract reads from the padded region
                    for read in bam.fetch(chrom, padded_start, padded_end):
                        # Skip unmapped and low-quality reads
                        if read.is_unmapped or read.mapping_quality < 10:
                            continue
                        
                        # Track unique reads (avoid duplicates from supplementary alignments)
                        if read.query_name in read_names:
                            continue
                        read_names.add(read.query_name)
                        
                        # Get read sequence and quality
                        sequence = read.query_sequence
                        quality = read.qual
                        
                        if not sequence:
                            continue
                        
                        # Write to FASTQ
                        fq.write(f"@{read.query_name}\n")
                        fq.write(f"{sequence}\n")
                        fq.write("+\n")
                        fq.write(f"{quality if quality else 'I' * len(sequence)}\n")
                        
                        # Update statistics
                        stats['reads_extracted'] += 1
                        read_lengths.append(len(sequence))
                        
                        if read.is_supplementary:
                            stats['supplementary_reads'] += 1
                        else:
                            stats['primary_reads'] += 1
                        
                        # Check for circular DNA signatures
                        if self._check_junction_spanning(read, start, end):
                            stats['reads_spanning_junction'] += 1
                        
                        if self._check_chimeric(read):
                            stats['chimeric_reads'] += 1
                    
                    # Calculate statistics
                    if read_lengths:
                        stats['avg_read_length'] = sum(read_lengths) / len(read_lengths)
                        stats['coverage_depth'] = sum(read_lengths) / stats['length']
                        stats['read_n50'] = self._calculate_n50(read_lengths)
            
            # Write detailed info file
            with open(info_file, 'w') as f:
                f.write(f"CircONTrack Region Information\n")
                f.write("=" * 50 + "\n")
                f.write(f"Region: {chrom}:{start}-{end}\n")
                f.write(f"Region ID: {region_id}\n")
                f.write(f"Region Type: {stats['region_type']}\n")
                f.write(f"Length: {stats['length']:,} bp\n")
                f.write(f"Padded region: {chrom}:{padded_start}-{padded_end}\n")
                f.write(f"Padding: {self.padding} bp\n")
                f.write("\nRead Statistics:\n")
                f.write(f"Reads extracted: {stats['reads_extracted']}\n")
                f.write(f"Primary reads: {stats['primary_reads']}\n")
                f.write(f"Supplementary reads: {stats['supplementary_reads']}\n")
                f.write(f"Junction-spanning reads: {stats['reads_spanning_junction']}\n")
                f.write(f"Chimeric reads: {stats['chimeric_reads']}\n")
                f.write(f"Average read length: {stats['avg_read_length']:.0f} bp\n")
                f.write(f"Read N50: {stats['read_n50']:,} bp\n")
                f.write(f"Estimated coverage: {stats['coverage_depth']:.1f}x\n")
                
                # Add assembly recommendations
                f.write("\nAssembly Recommendation:\n")
                if stats['coverage_depth'] < 10:
                    f.write("⚠️  Low coverage - assembly may be fragmented\n")
                elif stats['coverage_depth'] > 100:
                    f.write("✓ High coverage - good for assembly\n")
                else:
                    f.write("✓ Moderate coverage - suitable for assembly\n")
                
                if stats['reads_spanning_junction'] > 5:
                    f.write("✓ Junction-spanning reads detected - circular structure likely\n")
            
            return stats
            
        except Exception as e:
            print(f"Error processing region {region_id}: {str(e)}")
            stats['error'] = str(e)
            return stats
    
    def _check_junction_spanning(self, read, start, end):
        """Check if read might span a circular junction"""
        # Check for supplementary alignments (split alignments)
        if read.has_tag('SA'):
            sa_tag = read.get_tag('SA')
            # Parse SA tag to check if it maps back to same region
            for sa_entry in sa_tag.rstrip(';').split(';'):
                if sa_entry:
                    parts = sa_entry.split(',')
                    if len(parts) >= 2:
                        sa_pos = int(parts[1])
                        # Check if SA maps near the region boundaries
                        if abs(sa_pos - start) < 1000 or abs(sa_pos - end) < 1000:
                            return True
        
        # Check for large soft clips (potential junction signatures)
        if read.cigartuples:
            first_op = read.cigartuples[0]
            last_op = read.cigartuples[-1]
            
            # Large soft clips might indicate junction
            if (first_op[0] == 4 and first_op[1] > 100) or \
               (last_op[0] == 4 and last_op[1] > 100):
                return True
        
        return False
    
    def _check_chimeric(self, read):
        """Check if read has chimeric alignment signature"""
        # Has supplementary alignment
        if read.has_tag('SA'):
            return True
        # Is supplementary alignment itself
        if read.is_supplementary:
            return True
        return False
    
    def _calculate_n50(self, lengths):
        """Calculate N50 of read lengths"""
        if not lengths:
            return 0
        sorted_lengths = sorted(lengths, reverse=True)
        total = sum(sorted_lengths)
        target = total / 2
        cumsum = 0
        for length in sorted_lengths:
            cumsum += length
            if cumsum >= target:
                return length
        return sorted_lengths[-1]
    
    def process_bed_file(self, bed_file, classified=False):
        """
        Process all regions from a BED file
        
        Parameters:
        -----------
        bed_file : str
            Path to BED file (CircONTrack output or classified output)
        classified : bool
            Whether the BED file is from circontrack-classify (has type column)
        """
        # Read BED file
        regions = []
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    region = {
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3],
                        'region_id': f"{parts[0]}_{parts[1]}_{parts[2]}_{parts[3]}"
                    }
                    
                    # If classified BED, extract region type
                    if classified and len(parts) >= 7:
                        region['type'] = parts[6]
                    else:
                        region['type'] = None
                    
                    regions.append(region)
        
        print(f"Processing {len(regions)} circular DNA regions...")
        
        # Process regions in parallel
        all_stats = []
        
        if self.threads > 1:
            with ThreadPoolExecutor(max_workers=self.threads) as executor:
                futures = {
                    executor.submit(
                        self.extract_reads_from_region,
                        region['chrom'],
                        region['start'],
                        region['end'],
                        region['region_id'],
                        region.get('type')
                    ): region for region in regions
                }
                
                # Progress bar
                for future in tqdm(as_completed(futures), total=len(futures), 
                                  desc="Extracting reads"):
                    stats = future.result()
                    all_stats.append(stats)
        else:
            # Single-threaded processing
            for region in tqdm(regions, desc="Extracting reads"):
                stats = self.extract_reads_from_region(
                    region['chrom'],
                    region['start'],
                    region['end'],
                    region['region_id'],
                    region.get('type')
                )
                all_stats.append(stats)
        
        # Write summary statistics
        self._write_summary(all_stats)
        
        return all_stats
    
    def _write_summary(self, all_stats):
        """Write comprehensive summary statistics"""
        summary_file = self.output_dir / "extraction_summary.txt"
        stats_file = self.output_dir / "extraction_stats.tsv"
        
        # Calculate totals
        total_reads = sum(s['reads_extracted'] for s in all_stats)
        total_regions = len(all_stats)
        successful_regions = sum(1 for s in all_stats if 'error' not in s)
        
        # Group by type if available
        by_type = defaultdict(list)
        for s in all_stats:
            by_type[s.get('region_type', 'unknown')].append(s)
        
        # Write summary
        with open(summary_file, 'w') as f:
            f.write("CircONTrack Assembly Preparation Summary\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Total regions processed: {total_regions}\n")
            f.write(f"Successful extractions: {successful_regions}\n")
            f.write(f"Failed extractions: {total_regions - successful_regions}\n")
            f.write(f"Total reads extracted: {total_reads:,}\n")
            f.write(f"Output directory: {self.output_dir}\n")
            
            # Summary by type
            if len(by_type) > 1:
                f.write("\nRegions by type:\n")
                for region_type, stats_list in sorted(by_type.items()):
                    f.write(f"  {region_type}: {len(stats_list)} regions, "
                           f"{sum(s['reads_extracted'] for s in stats_list):,} reads\n")
            
            # Top regions by coverage
            f.write("\nTop 10 regions by coverage:\n")
            sorted_stats = sorted(all_stats, 
                                 key=lambda x: x.get('coverage_depth', 0), 
                                 reverse=True)
            for s in sorted_stats[:10]:
                f.write(f"  {s['region_id']}: {s['coverage_depth']:.1f}x coverage, "
                       f"{s['reads_extracted']:,} reads\n")
            
            # Regions with junction-spanning reads
            junction_regions = [s for s in all_stats 
                               if s.get('reads_spanning_junction', 0) > 0]
            if junction_regions:
                f.write(f"\nRegions with junction-spanning reads: {len(junction_regions)}\n")
                for s in sorted(junction_regions, 
                               key=lambda x: x['reads_spanning_junction'], 
                               reverse=True)[:5]:
                    f.write(f"  {s['region_id']}: {s['reads_spanning_junction']} junction reads\n")
        
        # Write detailed stats table
        if HAS_PANDAS:
            df = pd.DataFrame(all_stats)
            df.to_csv(stats_file, sep='\t', index=False)
        else:
            # Write TSV manually
            with open(stats_file, 'w') as f:
                # Header
                if all_stats:
                    keys = all_stats[0].keys()
                    f.write('\t'.join(str(k) for k in keys) + '\n')
                    # Data
                    for stats in all_stats:
                        f.write('\t'.join(str(stats.get(k, '')) for k in keys) + '\n')
        
        print(f"\nSummary written to: {summary_file}")
        print(f"Detailed stats written to: {stats_file}")
    
    def create_assembly_script(self, min_reads=50, min_coverage=10, 
                              assembler='flye', use_slurm=False):
        """
        Create assembly script for extracted reads
        
        Parameters:
        -----------
        min_reads : int
            Minimum number of reads to attempt assembly
        min_coverage : float
            Minimum coverage depth to attempt assembly
        assembler : str
            Assembler to use ('flye', 'canu', 'unicycler')
        use_slurm : bool
            Create SLURM/Swarm submission script
        """
        # Get all FASTQ files and filter by thresholds
        fastq_files = list(self.reads_dir.glob("*.fastq"))
        valid_jobs = []
        
        # Read statistics to determine which regions to assemble
        stats_file = self.output_dir / "extraction_stats.tsv"
        stats_dict = {}
        
        if stats_file.exists():
            with open(stats_file) as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == len(header):
                        stats = dict(zip(header, parts))
                        stats_dict[stats['region_id']] = stats
        
        # Filter regions by thresholds
        for fq in fastq_files:
            region_id = fq.stem
            
            if region_id in stats_dict:
                stats = stats_dict[region_id]
                reads = int(stats.get('reads_extracted', 0))
                coverage = float(stats.get('coverage_depth', 0))
                
                if reads >= min_reads and coverage >= min_coverage:
                    valid_jobs.append({
                        'region_id': region_id,
                        'fastq_file': fq,
                        'reads': reads,
                        'coverage': coverage,
                        'region_type': stats.get('region_type', 'unknown')
                    })
        
        if assembler == 'flye':
            self._create_flye_script(valid_jobs, use_slurm)
        elif assembler == 'canu':
            self._create_canu_script(valid_jobs, use_slurm)
        else:
            print(f"Assembler '{assembler}' not yet supported")
            return
        
        # Write job summary
        job_summary = self.output_dir / "assembly_jobs.txt"
        with open(job_summary, 'w') as f:
            f.write("CircONTrack Assembly Jobs\n")
            f.write("=" * 40 + "\n")
            f.write(f"Assembler: {assembler}\n")
            f.write(f"Total FASTQ files: {len(fastq_files)}\n")
            f.write(f"Jobs meeting criteria: {len(valid_jobs)}\n")
            f.write(f"Minimum reads: {min_reads}\n")
            f.write(f"Minimum coverage: {min_coverage}x\n\n")
            
            if valid_jobs:
                f.write("Regions to assemble:\n")
                for job in sorted(valid_jobs, key=lambda x: x['coverage'], reverse=True):
                    f.write(f"  {job['region_id']}: {job['reads']:,} reads, "
                           f"{job['coverage']:.1f}x coverage, "
                           f"type={job['region_type']}\n")
            else:
                f.write("No regions meet the assembly criteria.\n")
        
        print(f"\nAssembly jobs summary: {job_summary}")
        print(f"Regions to assemble: {len(valid_jobs)}/{len(fastq_files)}")
        
        return valid_jobs
    
    def _create_flye_script(self, valid_jobs, use_slurm):
        """Create Flye assembly script"""
        if use_slurm:
            # Create SLURM/Swarm file
            swarm_file = self.output_dir / "flye_assembly.swarm"
            submit_script = self.output_dir / "submit_flye.sh"
            
            with open(swarm_file, 'w') as f:
                f.write("# Flye assembly for CircONTrack regions\n")
                for job in valid_jobs:
                    region_id = job['region_id']
                    fastq = job['fastq_file']
                    out_dir = self.output_dir / f"assembly_{region_id}"
                    
                    cmd = (f"flye --nano-raw {fastq} "
                          f"--out-dir {out_dir} "
                          f"--threads $SLURM_CPUS_PER_TASK "
                          f"--meta --min-overlap 1000 --iterations 2")
                    f.write(f"{cmd}\n")
            
            with open(submit_script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(f"swarm -f {swarm_file} -t {self.threads} -g 16 "
                       f"--time 4:00:00 --module flye\n")
            
            os.chmod(submit_script, 0o755)
            print(f"SLURM submission script: {submit_script}")
        
        else:
            # Create regular bash script
            script_file = self.output_dir / "run_flye_assembly.sh"
            
            with open(script_file, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("# Flye assembly script for CircONTrack regions\n\n")
                
                for job in valid_jobs:
                    region_id = job['region_id']
                    fastq = job['fastq_file']
                    out_dir = self.output_dir / f"assembly_{region_id}"
                    
                    f.write(f"echo 'Assembling {region_id}...'\n")
                    f.write(f"flye --nano-raw {fastq} "
                           f"--out-dir {out_dir} "
                           f"--threads {self.threads} "
                           f"--meta --min-overlap 1000 --iterations 2\n\n")
            
            os.chmod(script_file, 0o755)
            print(f"Assembly script: {script_file}")
    
    def _create_canu_script(self, valid_jobs, use_slurm):
        """Create Canu assembly script"""
        script_file = self.output_dir / "run_canu_assembly.sh"
        
        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Canu assembly script for CircONTrack regions\n\n")
            
            for job in valid_jobs:
                region_id = job['region_id']
                fastq = job['fastq_file']
                out_dir = self.output_dir / f"assembly_{region_id}"
                genome_size = job.get('length', 10000)  # Estimate
                
                f.write(f"echo 'Assembling {region_id}...'\n")
                f.write(f"canu -p {region_id} -d {out_dir} "
                       f"genomeSize={genome_size} "
                       f"-nanopore {fastq}\n\n")
        
        os.chmod(script_file, 0o755)
        print(f"Assembly script: {script_file}")


def main():
    """Main entry point for circontrack-assemble"""
    parser = argparse.ArgumentParser(
        description="CircONTrack Assemble - Extract reads and prepare for assembly",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This module extracts reads from circular DNA regions detected by CircONTrack
and prepares them for assembly with tools like Flye or Canu.

Examples:
  # Basic usage
  circontrack-assemble circdna.bed sample.bam -o assembly_prep
  
  # Use classified BED with region types
  circontrack-assemble classified.bed sample.bam --classified -o assembly_prep
  
  # Custom thresholds and assembler
  circontrack-assemble circdna.bed sample.bam \\
    --min-reads 100 --min-coverage 20 \\
    --assembler flye --slurm \\
    -o assembly_prep
        """
    )
    
    parser.add_argument('bed_file', 
                       help='BED file with circular DNA regions (from CircONTrack)')
    parser.add_argument('bam_file', 
                       help='BAM file with aligned reads')
    
    parser.add_argument('-o', '--output-dir', default='assembly_prep',
                       help='Output directory (default: assembly_prep)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads (default: 4)')
    parser.add_argument('-p', '--padding', type=int, default=1000,
                       help='Padding around regions in bp (default: 1000)')
    
    parser.add_argument('--classified', action='store_true',
                       help='BED file is from circontrack-classify (has type column)')
    
    # Assembly parameters
    parser.add_argument('--min-reads', type=int, default=50,
                       help='Minimum reads for assembly (default: 50)')
    parser.add_argument('--min-coverage', type=float, default=10,
                       help='Minimum coverage for assembly (default: 10)')
    parser.add_argument('--assembler', choices=['flye', 'canu'], default='flye',
                       help='Assembler to use (default: flye)')
    parser.add_argument('--slurm', action='store_true',
                       help='Create SLURM/Swarm submission script')
    
    args = parser.parse_args()
    
    # Validate input files
    for filepath, name in [(args.bed_file, 'BED file'),
                           (args.bam_file, 'BAM file')]:
        if not os.path.exists(filepath):
            print(f"Error: {name} not found: {filepath}")
            sys.exit(1)
    
    print("CircONTrack Assemble Module")
    print("=" * 40)
    
    # Create assembler instance
    assembler = CircularDNAAssembler(
        args.bam_file,
        args.output_dir,
        args.threads,
        args.padding
    )
    
    # Process regions
    stats = assembler.process_bed_file(args.bed_file, args.classified)
    
    # Create assembly scripts
    if stats:
        jobs = assembler.create_assembly_script(
            args.min_reads,
            args.min_coverage,
            args.assembler,
            args.slurm
        )
        
        print("\n" + "=" * 40)
        print("Extraction complete!")
        print(f"Reads extracted to: {args.output_dir}/reads/")
        print(f"Region info in: {args.output_dir}/info/")
        
        if jobs:
            if args.slurm:
                print(f"\nSubmit assembly jobs with:")
                print(f"  bash {args.output_dir}/submit_{args.assembler}.sh")
            else:
                print(f"\nRun assembly with:")
                print(f"  bash {args.output_dir}/run_{args.assembler}_assembly.sh")
        else:
            print("\nNo regions meet assembly criteria.")
            print(f"Adjust --min-reads or --min-coverage if needed.")


if __name__ == "__main__":
    main()