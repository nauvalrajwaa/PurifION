#!/usr/bin/env python3
"""
downsampler.py — Smart Uniform Downsampler for Circular Topologies
===================================================================
Performs coverage-normalized subsampling of ONT reads for plant mitogenome
assembly, prioritising:
  1. Longest reads
  2. Highest Phred Q-score (average base quality)

Multi-chromosome support
------------------------
Plant mitogenomes commonly consist of multiple chromosomes (2–7 is typical).
The --reference FASTA may contain any number of sequences.  Bin-based coverage
normalisation, spanning-read detection, and coverage diagnostics all operate
per contig so that every chromosome independently reaches the target coverage,
regardless of differences in chromosome size or sequencing depth.

Algorithm
---------
1. Map all reads to the mitochondrial reference (minimap2).
2. Divide each contig into non-overlapping bins (default: 500 bp).
3. For each bin, collect all reads that *overlap* the bin.
4. Rank reads by (length DESC, mean_q DESC) and greedily select reads until
   the bin reaches target_coverage × bin_size bases of coverage.  Selected
   reads are added to a global "keep" set (reads covering multiple bins count
   for all of them).
5. Spanning reads — reads whose alignment covers ≥ spanning_fraction of their
   own contig's length — are always kept to resolve circular topology of each
   chromosome independently.
6. Reads in the keep set are written to the output file.

Usage
-----
    python downsampler.py \\
        --reads clean_reads.fastq \\
        --reference mito_ref.fasta \\
        --output downsampled.fastq \\
        [--target-coverage 100] \\
        [--bin-size 500] \\
        [--threads 8] \\
        [--spanning-fraction 0.80] \\
        [--stats stats.json]
"""

import argparse
import json
import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List, Set, Tuple

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Mapping
# ---------------------------------------------------------------------------

def run_minimap2(
    reference: str,
    reads: str,
    bam_out: str,
    threads: int = 4,
    preset: str = "map-ont",
) -> None:
    """Align reads and write sorted, indexed BAM."""
    sam_tmp = bam_out + ".tmp.sam"
    cmd_map = [
        "minimap2",
        "-a",
        "-x", preset,
        "-t", str(threads),
        "--secondary=no",
        reference,
        reads,
    ]
    log.info("minimap2: %s", " ".join(cmd_map))
    with open(sam_tmp, "w") as fh:
        subprocess.run(cmd_map, stdout=fh, check=True)

    subprocess.run(
        ["samtools", "sort", "-o", bam_out, "-@", str(threads), sam_tmp], check=True
    )
    os.remove(sam_tmp)
    subprocess.run(["samtools", "index", bam_out], check=True)


# ---------------------------------------------------------------------------
# Read info extraction
# ---------------------------------------------------------------------------

class ReadInfo:
    """Lightweight container for per-read statistics."""

    __slots__ = ("name", "length", "mean_q", "ref_name", "ref_start", "ref_end")

    def __init__(
        self,
        name: str,
        length: int,
        mean_q: float,
        ref_name: str,
        ref_start: int,
        ref_end: int,
    ) -> None:
        self.name = name
        self.length = length
        self.mean_q = mean_q
        self.ref_name = ref_name
        self.ref_start = ref_start
        self.ref_end = ref_end


def load_read_info(bam_path: str, min_mapq: int = 10) -> List[ReadInfo]:
    """
    Parse BAM and return list of ReadInfo for primary alignments with
    MAPQ ≥ min_mapq.

    Supports multi-contig references (multi-chromosome plant mitogenomes).
    Each read's ref_name, ref_start, ref_end are preserved so downstream
    functions can operate per-contig.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    reads: List[ReadInfo] = []

    for aln in bam.fetch():
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        query_len = aln.infer_read_length() or aln.query_length or 1
        # Mean base quality from the aligned portion; fall back to 0
        quals = aln.query_qualities
        mean_q = float(sum(quals)) / len(quals) if quals is not None and len(quals) > 0 else 0.0

        reads.append(
            ReadInfo(
                name=aln.query_name,
                length=query_len,
                mean_q=mean_q,
                ref_name=aln.reference_name,
                ref_start=aln.reference_start,
                ref_end=aln.reference_end,
            )
        )

    bam.close()
    return reads


# ---------------------------------------------------------------------------
# Spanning read detection
# ---------------------------------------------------------------------------

def find_spanning_reads(
    reads: List[ReadInfo],
    ref_lengths: Dict[str, int],
    spanning_fraction: float = 0.80,
) -> Set[str]:
    """
    Return names of reads whose aligned length covers ≥ spanning_fraction of
    their own contig's length.

    These reads are always retained regardless of the bin budget because they
    are critical for resolving the circular topology of each chromosome.

    Multi-chromosome note: spanning fraction is evaluated independently per
    contig.  A read spanning 80% of a 50 kb chromosome is kept for that
    chromosome; it does not need to span 80% of the total genome length.
    This is the correct behaviour — each chromosome needs its own spanning
    reads to be assembled into a closed circle.
    """
    spanning: Set[str] = set()
    for r in reads:
        contig_len = ref_lengths.get(r.ref_name, 1)
        aligned_len = r.ref_end - r.ref_start
        if aligned_len / contig_len >= spanning_fraction:
            spanning.add(r.name)
    return spanning


# ---------------------------------------------------------------------------
# Bin-based greedy selection
# ---------------------------------------------------------------------------

def greedy_bin_selection(
    reads: List[ReadInfo],
    ref_lengths: Dict[str, int],
    target_coverage: int = 100,
    bin_size: int = 500,
) -> Set[str]:
    """
    Greedy bin-filling selection.

    For each bin across the reference, rank candidate reads by
    (length DESC, mean_q DESC) and add them to the keep set until
    the bin reaches target_coverage × bin_size bases of coverage.

    Returns the set of selected read names.
    """
    # Build bin → reads map
    # bin key: (ref_name, bin_index)
    bin_reads: Dict[Tuple[str, int], List[ReadInfo]] = defaultdict(list)

    for r in reads:
        contig_len = ref_lengths.get(r.ref_name, 1)
        n_bins = max(1, (contig_len + bin_size - 1) // bin_size)
        # Find all bins overlapping this alignment
        first_bin = r.ref_start // bin_size
        last_bin = min((r.ref_end - 1) // bin_size, n_bins - 1)
        for b in range(first_bin, last_bin + 1):
            bin_reads[(r.ref_name, b)].append(r)

    # Sort each bin's candidate list: longest first, then highest Q
    for key in bin_reads:
        bin_reads[key].sort(key=lambda r: (r.length, r.mean_q), reverse=True)

    keep: Set[str] = set()
    # Track accumulated coverage per bin (in bases)
    bin_coverage: Dict[Tuple[str, int], int] = defaultdict(int)
    # How many bases of coverage we want per bin
    coverage_budget = target_coverage * bin_size

    for (ref_name, bin_idx), candidates in sorted(bin_reads.items()):
        for r in candidates:
            if bin_coverage[(ref_name, bin_idx)] >= coverage_budget:
                break
            keep.add(r.name)
            # Credit coverage to ALL bins this read overlaps
            contig_len = ref_lengths.get(ref_name, 1)
            n_bins = max(1, (contig_len + bin_size - 1) // bin_size)
            first_b = r.ref_start // bin_size
            last_b = min((r.ref_end - 1) // bin_size, n_bins - 1)
            overlap = r.ref_end - r.ref_start  # approximate; good enough for budgeting
            per_bin_bases = overlap // max(last_b - first_b + 1, 1)
            for b in range(first_b, last_b + 1):
                bin_coverage[(ref_name, b)] += per_bin_bases

    return keep


# ---------------------------------------------------------------------------
# Read extraction
# ---------------------------------------------------------------------------

def _detect_format(path: str) -> str:
    with open(path) as fh:
        first = fh.read(1)
    return "fastq" if first == "@" else "fasta"


def extract_reads(
    reads_path: str,
    keep_names: Set[str],
    output_path: str,
) -> int:
    """Write reads in keep_names to output_path. Returns count written."""
    fmt = _detect_format(reads_path)
    written = 0
    with pysam.FastxFile(reads_path) as fin, open(output_path, "w") as fout:
        for entry in fin:
            if entry.name in keep_names:
                if fmt == "fastq":
                    fout.write(
                        f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}\n"
                    )
                else:
                    fout.write(f">{entry.name}\n{entry.sequence}\n")
                written += 1
    return written


# ---------------------------------------------------------------------------
# Coverage diagnostics
# ---------------------------------------------------------------------------

def compute_achieved_coverage(
    keep_names: Set[str],
    reads: List[ReadInfo],
    ref_lengths: Dict[str, int],
) -> Dict[str, float]:
    """Return dict of {contig_name: achieved_mean_coverage} for kept reads."""
    cov_sum: Dict[str, int] = defaultdict(int)
    for r in reads:
        if r.name in keep_names:
            cov_sum[r.ref_name] += r.ref_end - r.ref_start
    return {
        ref: cov_sum.get(ref, 0) / length
        for ref, length in ref_lengths.items()
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Smart Uniform Downsampler for plant mitogenome circular topology.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--reads", required=True, help="Input reads (FASTQ or FASTA)")
    p.add_argument("--reference", required=True, help="Mitochondrial reference FASTA")
    p.add_argument("--output", required=True, help="Output downsampled reads")
    p.add_argument("--target-coverage", type=int, default=100, dest="target_coverage")
    p.add_argument("--bin-size", type=int, default=500, dest="bin_size",
                   help="Reference bin size (bp) for coverage normalisation")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-mapq", type=int, default=10, dest="min_mapq")
    p.add_argument(
        "--spanning-fraction",
        type=float,
        default=0.80,
        dest="spanning_fraction",
        help="Reads covering ≥ this fraction of the genome are always kept",
    )
    p.add_argument("--stats", default=None, help="Optional path for JSON statistics")
    p.add_argument(
        "--keep-bam",
        action="store_true",
        dest="keep_bam",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = os.path.join(tmpdir, "aligned.bam")

        # 1. Map
        run_minimap2(
            reference=args.reference,
            reads=args.reads,
            bam_out=bam_path,
            threads=args.threads,
        )

        if args.keep_bam:
            import shutil
            kept = args.output + ".downsampler.bam"
            shutil.copy(bam_path, kept)
            shutil.copy(bam_path + ".bai", kept + ".bai")
            log.info("BAM kept at: %s", kept)

        # 2. Load read info
        log.info("Loading alignment info …")
        reads = load_read_info(bam_path, min_mapq=args.min_mapq)
        log.info("Loaded %d primary alignments", len(reads))

        bam = pysam.AlignmentFile(bam_path, "rb")
        ref_lengths: Dict[str, int] = dict(zip(bam.references, bam.lengths))
        bam.close()

        # 3. Always-keep spanning reads
        spanning = find_spanning_reads(reads, ref_lengths, args.spanning_fraction)
        log.info("Spanning reads (always kept): %d", len(spanning))

        # 4. Greedy bin selection
        log.info("Running greedy bin selection (target=%dx, bin=%dbp) …",
                 args.target_coverage, args.bin_size)
        selected = greedy_bin_selection(
            reads,
            ref_lengths,
            target_coverage=args.target_coverage,
            bin_size=args.bin_size,
        )
        keep = selected | spanning
        log.info("Total reads selected: %d (bin=%d + spanning_only=%d)",
                 len(keep), len(selected), len(spanning - selected))

        # 5. Extract
        written = extract_reads(args.reads, keep, args.output)
        log.info("Wrote %d reads to %s", written, args.output)

        # 6. Coverage diagnostics
        achieved = compute_achieved_coverage(keep, reads, ref_lengths)
        for ref, cov in achieved.items():
            log.info("Achieved coverage on %s: %.1fx", ref, cov)

        # 7. Stats
        stats = {
            "input_alignments": len(reads),
            "spanning_reads": len(spanning),
            "selected_reads": len(keep),
            "output_reads": written,
            "target_coverage": args.target_coverage,
            "achieved_coverage": achieved,
        }
        if args.stats:
            with open(args.stats, "w") as fh:
                json.dump(stats, fh, indent=2)
            log.info("Statistics written to %s", args.stats)


if __name__ == "__main__":
    main()
