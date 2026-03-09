#!/usr/bin/env python3
"""
numt_scrubber.py — Context & Depth-Aware NUMT Scrubber
=======================================================
Filters ONT long reads to remove NUMTs (Nuclear Mitochondrial DNA segments)
and chimeric reads (part mitochondrial, part nuclear) before assembly.

Multi-chromosome support
------------------------
Plant mitogenomes often consist of multiple chromosomes (2–7 is common).
The --reference FASTA may contain any number of sequences.  Coverage baselines
and thresholds are computed **per contig** so that a chromosome with 200x
coverage and one with 50x coverage each get their own independent threshold,
preventing reads on the low-coverage chromosome from being unfairly discarded.

Algorithm
---------
1. Map all reads to the mitochondrial reference (minimap2, map-ont preset).
2. Compute per-base coverage across every contig using a sliding window.
3. Estimate a per-contig baseline as the median of that contig's window values.
4. For each primary alignment:
   a. Coverage context: check whether the read maps to a region with
      coverage ≥ contig_baseline × coverage_threshold_ratio.  Low-coverage
      regions suggest NUMT insertion sites; threshold is contig-specific.
   b. Alignment quality: check that mapped fraction ≥ min_mapped_fraction
      (default: 0.80) so chimeric reads with large soft-clipped nuclear
      tails are removed.
   c. Flanking context: reads whose soft-clip exceeds max_softclip_fraction
      (default: 0.20 of read length) are flagged as chimeric/NUMT.
5. Reads passing all filters are written to the output FASTQ/FASTA.

Usage
-----
    python numt_scrubber.py \\
        --reads reads.fastq \\
        --reference mito_ref.fasta \\
        --output clean_reads.fastq \\
        [--threads 8] \\
        [--min-mapq 20] \\
        [--min-mapped-fraction 0.80] \\
        [--max-softclip-fraction 0.20] \\
        [--coverage-window 500] \\
        [--coverage-threshold-ratio 0.20] \\
        [--stats stats.json]
"""

import argparse
import json
import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from statistics import median
from typing import Dict, List, Optional, Set, Tuple

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def run_minimap2(
    reference: str,
    reads: str,
    bam_out: str,
    threads: int = 4,
    preset: str = "map-ont",
) -> None:
    """Align reads to reference with minimap2 and produce a sorted BAM."""
    sam_tmp = bam_out + ".unsorted.sam"
    cmd_map = [
        "minimap2",
        "-a",          # SAM output
        "-x", preset,
        "-t", str(threads),
        "--secondary=no",
        reference,
        reads,
    ]
    cmd_sort = [
        "samtools", "sort",
        "-o", bam_out,
        "-@", str(threads),
        sam_tmp,
    ]
    cmd_index = ["samtools", "index", bam_out]

    log.info("Running minimap2: %s", " ".join(cmd_map))
    with open(sam_tmp, "w") as fh:
        subprocess.run(cmd_map, stdout=fh, check=True)

    log.info("Sorting BAM …")
    subprocess.run(cmd_sort, check=True)
    os.remove(sam_tmp)

    log.info("Indexing BAM …")
    subprocess.run(cmd_index, check=True)


# ---------------------------------------------------------------------------
# Coverage analysis
# ---------------------------------------------------------------------------

def compute_sliding_window_coverage(
    bam_path: str,
    window: int = 500,
    step: int = 100,
) -> Dict[str, List[float]]:
    """
    Return per-contig list of mean coverage values in sliding windows.
    Keys are contig names; values are lists of (window_mean_coverage).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    contig_cov: Dict[str, List[float]] = {}

    for ref in bam.references:
        length = bam.get_reference_length(ref)
        # pysam count_coverage returns 4 arrays (A, C, G, T) per position
        cov_arrays = bam.count_coverage(ref, quality_threshold=0)
        total_cov = [
            sum(cov_arrays[b][i] for b in range(4))
            for i in range(length)
        ]
        windows: List[float] = []
        for start in range(0, length - window + 1, step):
            end = start + window
            windows.append(sum(total_cov[start:end]) / window)
        if not windows and length > 0:
            # contig shorter than window — single window
            windows.append(sum(total_cov) / length)
        contig_cov[ref] = windows

    bam.close()
    return contig_cov


def estimate_baselines(
    cov_data: Dict[str, List[float]],
) -> Dict[str, float]:
    """
    Return a per-contig dict of {contig_name: median_window_coverage}.

    Each contig in a multi-chromosome mitogenome gets its own independent
    baseline so that coverage differences between chromosomes do not cause
    unfair filtering of reads on lower-coverage contigs.
    """
    return {
        contig: (median(windows) if windows else 0.0)
        for contig, windows in cov_data.items()
    }


# Backward-compat alias — kept so existing callers that pass a flat list still work.
def estimate_baseline(window_coverages: List[float]) -> float:
    """Return median coverage over all windows as a single scalar baseline.

    Deprecated for multi-contig use; prefer estimate_baselines().
    """
    if not window_coverages:
        return 0.0
    return median(window_coverages)


# ---------------------------------------------------------------------------
# Per-read classification
# ---------------------------------------------------------------------------

def classify_reads(
    bam_path: str,
    per_contig_baselines: Dict[str, float],
    coverage_threshold_ratio: float = 0.20,
    min_mapq: int = 20,
    min_mapped_fraction: float = 0.80,
    max_softclip_fraction: float = 0.20,
) -> Tuple[Set[str], Dict]:
    """
    Return (set_of_clean_read_names, stats_dict).

    A read is CLEAN when ALL of the following hold:
    - MAPQ ≥ min_mapq
    - mapped fraction of query ≥ min_mapped_fraction
    - soft-clip fraction ≤ max_softclip_fraction
    - maps to a region whose coverage ≥ contig_baseline × coverage_threshold_ratio

    per_contig_baselines: dict mapping contig name → median coverage baseline.
    Each contig is evaluated against its own baseline so multi-chromosome
    references (e.g. plant mitogenomes with 2–7 chromosomes) are handled
    correctly regardless of coverage differences between chromosomes.
    """
    stats = {
        "total_alignments": 0,
        "failed_mapq": 0,
        "failed_mapped_fraction": 0,
        "failed_softclip": 0,
        "failed_low_coverage": 0,
        "clean_reads": 0,
        "numt_or_chimeric_reads": 0,
        "per_contig_baselines": {k: round(v, 2) for k, v in per_contig_baselines.items()},
    }

    # Pre-compute per-contig thresholds
    per_contig_threshold: Dict[str, float] = {
        contig: baseline * coverage_threshold_ratio
        for contig, baseline in per_contig_baselines.items()
    }
    for contig, thr in per_contig_threshold.items():
        log.info(
            "  Contig %-30s  baseline=%.1f  threshold=%.1f (ratio=%.2f)",
            contig,
            per_contig_baselines[contig],
            thr,
            coverage_threshold_ratio,
        )

    bam = pysam.AlignmentFile(bam_path, "rb")
    # Pre-compute per-position total coverage for each contig
    per_contig_cov: Dict[str, List[int]] = {}
    for ref in bam.references:
        length = bam.get_reference_length(ref)
        cov_arrays = bam.count_coverage(ref, quality_threshold=0)
        per_contig_cov[ref] = [
            sum(cov_arrays[b][i] for b in range(4)) for i in range(length)
        ]

    clean: Set[str] = set()
    dirty: Set[str] = set()

    for aln in bam.fetch():
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue

        stats["total_alignments"] += 1
        read_name = aln.query_name
        query_len = aln.infer_read_length() or aln.query_length or 1

        # --- Filter 1: mapping quality ---
        if aln.mapping_quality < min_mapq:
            stats["failed_mapq"] += 1
            dirty.add(read_name)
            continue

        # --- Filter 2: mapped fraction ---
        mapped_bases = aln.query_alignment_length  # excludes soft/hard clips
        mapped_frac = mapped_bases / query_len
        if mapped_frac < min_mapped_fraction:
            stats["failed_mapped_fraction"] += 1
            dirty.add(read_name)
            continue

        # --- Filter 3: soft-clip fraction (chimeric indicator) ---
        cigar = aln.cigartuples or []
        softclip = sum(length for op, length in cigar if op == 4)  # op 4 = SOFT_CLIP
        softclip_frac = softclip / query_len
        if softclip_frac > max_softclip_fraction:
            stats["failed_softclip"] += 1
            dirty.add(read_name)
            continue

        # --- Filter 4: local coverage at alignment site ---
        # Threshold is looked up per-contig so multi-chromosome references
        # (different chromosomes may have very different coverage levels) are
        # handled correctly.
        ref_name = aln.reference_name
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        cov_slice = per_contig_cov.get(ref_name, [])
        if cov_slice:
            region_cov = sum(cov_slice[ref_start:ref_end]) / max(ref_end - ref_start, 1)
        else:
            region_cov = 0.0

        contig_threshold = per_contig_threshold.get(ref_name, 0.0)
        if region_cov < contig_threshold:
            stats["failed_low_coverage"] += 1
            dirty.add(read_name)
            continue

        clean.add(read_name)

    bam.close()

    # Reads that failed ANY filter are dirty; reads that never mapped are not
    # in clean set either — we keep only confirmed-clean reads.
    stats["clean_reads"] = len(clean)
    stats["numt_or_chimeric_reads"] = len(dirty)
    return clean, stats


# ---------------------------------------------------------------------------
# Read extraction
# ---------------------------------------------------------------------------

def extract_clean_reads(
    reads_path: str,
    clean_names: Set[str],
    output_path: str,
) -> int:
    """
    Write reads whose name is in clean_names to output_path.
    Handles FASTQ and FASTA (auto-detected). Returns count written.
    """
    fmt = _detect_format(reads_path)
    written = 0
    with pysam.FastxFile(reads_path) as fin, open(output_path, "w") as fout:
        for entry in fin:
            if entry.name in clean_names:
                if fmt == "fastq":
                    fout.write(
                        f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}\n"
                    )
                else:
                    fout.write(f">{entry.name}\n{entry.sequence}\n")
                written += 1
    return written


def _detect_format(path: str) -> str:
    with open(path) as fh:
        first = fh.read(1)
    return "fastq" if first == "@" else "fasta"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Context & Depth-Aware NUMT Scrubber for plant mitogenome reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--reads", required=True, help="Input reads (FASTQ or FASTA)")
    p.add_argument("--reference", required=True, help="Mitochondrial reference FASTA")
    p.add_argument("--output", required=True, help="Output clean reads (FASTQ/FASTA)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-mapq", type=int, default=20, dest="min_mapq")
    p.add_argument(
        "--min-mapped-fraction",
        type=float,
        default=0.80,
        dest="min_mapped_fraction",
        help="Minimum fraction of read that must align (0–1)",
    )
    p.add_argument(
        "--max-softclip-fraction",
        type=float,
        default=0.20,
        dest="max_softclip_fraction",
        help="Maximum soft-clip fraction allowed (0–1)",
    )
    p.add_argument(
        "--coverage-window",
        type=int,
        default=500,
        dest="coverage_window",
        help="Sliding window size (bp) for coverage estimation",
    )
    p.add_argument(
        "--coverage-threshold-ratio",
        type=float,
        default=0.20,
        dest="coverage_threshold_ratio",
        help="Fraction of baseline coverage below which a region is flagged as NUMT",
    )
    p.add_argument(
        "--stats",
        default=None,
        help="Optional path to write JSON statistics",
    )
    p.add_argument(
        "--keep-bam",
        action="store_true",
        dest="keep_bam",
        help="Keep intermediate BAM file instead of deleting it",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = os.path.join(tmpdir, "aligned.bam")

        # 1. Map reads
        run_minimap2(
            reference=args.reference,
            reads=args.reads,
            bam_out=bam_path,
            threads=args.threads,
        )

        if args.keep_bam:
            import shutil
            kept = args.output + ".numt_scrubber.bam"
            shutil.copy(bam_path, kept)
            shutil.copy(bam_path + ".bai", kept + ".bai")
            log.info("BAM kept at: %s", kept)

        # 2. Compute per-contig coverage baselines
        # Each mitochondrial chromosome gets its own baseline so multi-chromosome
        # plant mitogenomes (e.g. 2–7 chromosomes) are handled correctly.
        log.info("Computing sliding-window coverage …")
        cov_data = compute_sliding_window_coverage(
            bam_path, window=args.coverage_window
        )
        per_contig_baselines = estimate_baselines(cov_data)
        for contig, bl in per_contig_baselines.items():
            log.info("  Contig %-30s  baseline=%.1fx", contig, bl)

        # 3. Classify reads
        log.info("Classifying reads …")
        clean_names, stats = classify_reads(
            bam_path=bam_path,
            per_contig_baselines=per_contig_baselines,
            coverage_threshold_ratio=args.coverage_threshold_ratio,
            min_mapq=args.min_mapq,
            min_mapped_fraction=args.min_mapped_fraction,
            max_softclip_fraction=args.max_softclip_fraction,
        )

        # 4. Extract clean reads
        log.info("Extracting %d clean reads …", len(clean_names))
        written = extract_clean_reads(args.reads, clean_names, args.output)
        log.info("Wrote %d reads to %s", written, args.output)

        stats["output_reads"] = written

        if args.stats:
            with open(args.stats, "w") as fh:
                json.dump(stats, fh, indent=2)
            log.info("Statistics written to %s", args.stats)

        # Print summary
        total = stats["total_alignments"]
        log.info(
            "Summary: %d aligned | %d clean | %d NUMT/chimeric "
            "(mapq=%d mapped_frac=%d softclip=%d low_cov=%d)",
            total,
            stats["clean_reads"],
            stats["numt_or_chimeric_reads"],
            stats["failed_mapq"],
            stats["failed_mapped_fraction"],
            stats["failed_softclip"],
            stats["failed_low_coverage"],
        )


if __name__ == "__main__":
    main()
