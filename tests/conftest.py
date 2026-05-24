"""Deterministic synthetic genomics fixtures for CircONTrack tests."""

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pytest


Cigar = Sequence[Tuple[int, int]]


def require_pysam():
    """Skip pysam-dependent tests when the binary dependency is unavailable."""
    return pytest.importorskip("pysam")


def deterministic_sequence(length: int) -> str:
    """Return a reproducible DNA sequence of exactly ``length`` bases."""
    motif = "ACGT"
    return (motif * ((length // len(motif)) + 1))[:length]


def write_reference(tmp_path: Path, contigs: Optional[Dict[str, int]] = None) -> Path:
    """Create a tiny indexed FASTA reference matching synthetic BAM contigs."""
    pysam = require_pysam()
    contigs = contigs or {"chr1": 20000, "chr2": 12000}
    fasta_path = tmp_path / "reference.fa"
    with fasta_path.open("w") as handle:
        for name, length in contigs.items():
            handle.write(f">{name}\n")
            seq = deterministic_sequence(length)
            for i in range(0, length, 80):
                handle.write(seq[i:i + 80] + "\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


def build_header(contigs: Optional[Dict[str, int]] = None):
    """Build a minimal BAM header with deterministic contig ordering."""
    pysam = require_pysam()
    contigs = contigs or {"chr1": 20000, "chr2": 12000}
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs.items()],
    })


def query_length_from_cigar(cigar: Cigar) -> int:
    """Return query-consuming length for a pysam CIGAR tuple list."""
    return sum(length for op, length in cigar if op in {0, 1, 4, 7, 8})


def make_read(
    header,
    query_name: str,
    chrom: str,
    start: int,
    cigar: Cigar = ((0, 100),),
    mapq: int = 60,
    flag: int = 0,
    tags: Optional[Iterable[Tuple[str, object]]] = None,
):
    """Create one aligned read encoding a simple linear or clipped alignment."""
    pysam = require_pysam()
    read = pysam.AlignedSegment(header)
    read.query_name = query_name
    read.flag = flag
    read.reference_id = header.get_tid(chrom)
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigartuples = list(cigar)
    query_length = query_length_from_cigar(cigar)
    read.query_sequence = deterministic_sequence(query_length)
    read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    for tag, value in tags or []:
        read.set_tag(tag, value)
    return read


def make_split_junction_read(
    header,
    query_name: str,
    chrom: str = "chr1",
    primary_start: int = 100,
    supplementary_start: int = 250,
    segment_length: int = 150,
    mapq: int = 60,
):
    """Create a primary read with an SA tag supporting a back-to-back junction.

    The primary alignment is ``segment_length M`` followed by a soft clip. The
    SA tag maps the clipped segment back to the same contig at
    ``supplementary_start``. Coordinates in the SA tag are 1-based, as required
    by SAM, while pysam read coordinates remain 0-based.
    """
    primary_cigar = ((0, segment_length), (4, segment_length))
    supplementary_cigar = f"{segment_length}S{segment_length}M"
    sa_tag = f"{chrom},{supplementary_start + 1},+,{supplementary_cigar},{mapq},0;"
    return make_read(
        header,
        query_name=query_name,
        chrom=chrom,
        start=primary_start,
        cigar=primary_cigar,
        mapq=mapq,
        tags=[("SA", sa_tag), ("NM", 0)],
    )


def write_bam(
    tmp_path: Path,
    header,
    reads: Sequence,
    name: str = "synthetic",
    index: bool = True,
) -> Path:
    """Write, coordinate-sort, and optionally index a minimal synthetic BAM."""
    pysam = require_pysam()
    raw_path = tmp_path / f"{name}.unsorted.bam"
    bam_path = tmp_path / f"{name}.bam"

    with pysam.AlignmentFile(str(raw_path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)

    pysam.sort("-o", str(bam_path), str(raw_path))
    if index:
        pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def default_contigs() -> Dict[str, int]:
    return {"chr1": 20000, "chr2": 12000}
