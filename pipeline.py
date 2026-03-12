#!/usr/bin/env python3
"""
pipeline.py — MitoAssembler End-to-End Pipeline
================================================
Single-command wrapper that orchestrates:
  1. Per-chromosome PAF filter  (minimap2 PAF → per-chromosome read sets)
  2. Per-chromosome Downsampler (downsampler.py, one run per chromosome)
  3. Per-chromosome Flye assembly
  4. Disentangle graph          (disentangle_graph.py — mandatory)
  5. HTML Report                (report.py)

Usage
-----
    python pipeline.py \\
        --reads raw_reads.fastq \\
        --reference mito_ref.fasta \\
        --outdir results/ \\
        [--threads 16] \\
        [--target-coverage 100] \\
        [--run-assembler] \\
        [--assembler flye] \\
        [--genome-size 500k]

Output layout (inside --outdir)
--------------------------------
    01_per_chrom/
        <chrom>/
            filtered.fastq
            downsampled.fastq
            downsample_stats.json
    02_assembly/           (only when --run-assembler)
        <chrom>/
            assembly.fasta
            assembly_graph.gfa
            disentangle.report.txt
            disentangle.terminal_overlaps.tsv
            ...
    report/
        report.html
    pipeline_params.json
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam

# Local modules (same directory)
sys.path.insert(0, os.path.dirname(__file__))

import downsampler
import report as report_mod

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Reference helpers
# ---------------------------------------------------------------------------

def _get_chrom_names(reference: str) -> List[str]:
    """Return ordered list of sequence names from a FASTA reference."""
    names: List[str] = []
    with open(reference) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


# ---------------------------------------------------------------------------
# Per-chromosome PAF filter
# ---------------------------------------------------------------------------

def _run_paf_filter_per_chrom(
    reads: str,
    reference: str,
    chrom_names: List[str],
    out_dir: Path,
    threads: int,
    len_a: int = 5000,
    id_a: float = 0.90,
    len_b: int = 3000,
    id_b: float = 0.98,
) -> Dict[str, str]:
    """
    Align all reads against the reference with minimap2 (PAF output).
    For each chromosome, collect reads that pass the dual-threshold filter
    AND whose primary alignment target is that chromosome.

    Returns a dict: {chrom_name: path_to_filtered_fastq}
    """
    log.info(
        "PAF filter thresholds: A(len>=%d id>=%.2f)  B(len>=%d id>=%.2f)",
        len_a, id_a, len_b, id_b,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        paf_path = os.path.join(tmpdir, "aligned.paf")
        mm2_cmd = [
            "minimap2", "-cx", "map-ont",
            "-t", str(threads),
            reference, reads,
        ]
        log.info("minimap2 PAF: %s", " ".join(mm2_cmd))
        with open(paf_path, "w") as paf_fh:
            subprocess.run(mm2_cmd, stdout=paf_fh, check=True)

        # chrom → set of read names that passed filter targeting that chrom
        chrom_reads: Dict[str, Set[str]] = {c: set() for c in chrom_names}
        # Track best (score) alignment per read to assign it to one chrom
        read_best_score: Dict[str, Tuple[float, str]] = {}  # name → (score, chrom)

        total_aln = 0
        with open(paf_path) as paf_fh:
            for line in paf_fh:
                cols = line.split("\t")
                if len(cols) < 12:
                    continue
                total_aln += 1
                read_name = cols[0]
                target_name = cols[5]
                aln_len = int(cols[10])
                matches = int(cols[9])
                if aln_len == 0:
                    continue
                identity = matches / aln_len

                passes = (aln_len >= len_a and identity >= id_a) or \
                         (aln_len >= len_b and identity >= id_b)
                if not passes:
                    continue
                if target_name not in chrom_reads:
                    continue

                score = identity * aln_len
                prev = read_best_score.get(read_name)
                if prev is None or score > prev[0]:
                    # Remove from old chrom if re-assigned
                    if prev is not None:
                        chrom_reads[prev[1]].discard(read_name)
                    read_best_score[read_name] = (score, target_name)
                    chrom_reads[target_name].add(read_name)

        log.info(
            "PAF filter: %d alignments, reads per chrom: %s",
            total_aln,
            {c: len(v) for c, v in chrom_reads.items()},
        )

    # Write per-chromosome FASTQ files
    result: Dict[str, str] = {}
    # Accumulate all reads once into memory to avoid re-scanning for each chrom
    read_seqs: Dict[str, Tuple[str, str]] = {}  # name → (seq, qual)
    with pysam.FastxFile(reads) as fin:
        for entry in fin:
            if entry.name is not None:
                read_seqs[entry.name] = (entry.sequence or "", entry.quality or "")

    for chrom in chrom_names:
        chrom_dir = out_dir / chrom
        chrom_dir.mkdir(parents=True, exist_ok=True)
        out_fastq = str(chrom_dir / "filtered.fastq")
        written = 0
        with open(out_fastq, "w") as fout:
            for rname in chrom_reads[chrom]:
                if rname in read_seqs:
                    seq, qual = read_seqs[rname]
                    fout.write(f"@{rname}\n{seq}\n+\n{qual}\n")
                    written += 1
        log.info("  Chrom %-30s → %d reads → %s", chrom, written, out_fastq)
        result[chrom] = out_fastq

    return result


# ---------------------------------------------------------------------------
# Per-chromosome downsampler
# ---------------------------------------------------------------------------

def _downsample_chrom(
    filtered_reads: str,
    reference: str,
    chrom_name: str,
    out_fastq: str,
    stats_path: str,
    threads: int,
    target_coverage: int,
    bin_size: int,
    spanning_fraction: float,
) -> Dict:
    """
    Run downsampler for a single chromosome.
    Aligns filtered reads (already chromosome-specific) to full reference,
    restricts to the target chromosome, and runs greedy bin selection.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = os.path.join(tmpdir, "ds_aligned.bam")
        downsampler.run_minimap2(
            reference=reference,
            reads=filtered_reads,
            bam_out=bam_path,
            threads=threads,
        )
        ds_reads = downsampler.load_read_info(bam_path, min_mapq=10)

        bam = pysam.AlignmentFile(bam_path, "rb")
        ref_lengths = dict(zip(bam.references, bam.lengths))
        bam.close()

        # Restrict ref_lengths to just this chromosome so coverage is computed correctly
        chrom_ref_lengths = {chrom_name: ref_lengths[chrom_name]} if chrom_name in ref_lengths else ref_lengths

        # Filter ds_reads to only reads on this chromosome
        chrom_reads = [r for r in ds_reads if r.ref_name == chrom_name]
        if not chrom_reads:
            log.warning("  No mapped reads on %s after downsampler alignment", chrom_name)
            open(out_fastq, "w").close()
            stats: Dict = {
                "chrom": chrom_name,
                "input_alignments": 0,
                "spanning_reads": 0,
                "selected_reads": 0,
                "output_reads": 0,
                "target_coverage": target_coverage,
                "achieved_coverage": {},
                "ref_lengths": chrom_ref_lengths,
            }
            with open(stats_path, "w") as fh:
                json.dump(stats, fh, indent=2)
            return stats

        spanning = downsampler.find_spanning_reads(chrom_reads, chrom_ref_lengths, spanning_fraction)
        selected = downsampler.greedy_bin_selection(
            chrom_reads, chrom_ref_lengths, target_coverage=target_coverage, bin_size=bin_size
        )
        keep = selected | spanning

        ds_written = downsampler.extract_reads(filtered_reads, keep, out_fastq)
        log.info(
            "  Chrom %-30s downsampler: %d → %d reads",
            chrom_name, len(chrom_reads), ds_written,
        )

        achieved = downsampler.compute_achieved_coverage(keep, chrom_reads, chrom_ref_lengths)

    stats = {
        "chrom": chrom_name,
        "input_alignments": len(chrom_reads),
        "spanning_reads": len(spanning),
        "selected_reads": len(keep),
        "output_reads": ds_written,
        "target_coverage": target_coverage,
        "achieved_coverage": achieved,
        "ref_lengths": chrom_ref_lengths,
    }
    with open(stats_path, "w") as fh:
        json.dump(stats, fh, indent=2)
    return stats


# ---------------------------------------------------------------------------
# Assembler helpers
# ---------------------------------------------------------------------------

def _run_flye(
    reads: str,
    outdir: str,
    threads: int,
    genome_size: str,
    read_type: str = "--nano-hq",
) -> str:
    """Run Flye and return path to final assembly FASTA."""
    cmd = [
        "flye",
        read_type,
        reads,
        "--out-dir", outdir,
        "--threads", str(threads),
        "--genome-size", genome_size,
    ]
    log.info("Flye command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    assembly = os.path.join(outdir, "assembly.fasta")
    if not os.path.exists(assembly):
        raise FileNotFoundError(f"Flye finished but assembly.fasta not found in {outdir}")
    return assembly


ASSEMBLERS = {
    "flye": _run_flye,
}


def _run_disentangle(
    gfa_path: str,
    out_prefix: str,
    window: int,
    min_id: float,
    min_ovlp: int,
    top_out: int,
    expected_size: int,
    size_tolerance: float,
    max_depth: int,
    max_cycles: int,
    max_linear: int,
    max_cycles_any: int,
    fasta_top: int,
) -> Dict:
    """Run disentangle_graph.py on a GFA file (mandatory after assembly)."""
    disentangle_script = os.path.join(os.path.dirname(__file__), "disentangle_graph.py")
    if not os.path.exists(disentangle_script):
        raise FileNotFoundError(f"disentangle_graph.py not found at {disentangle_script}")

    cmd = [
        sys.executable,
        disentangle_script,
        gfa_path,
        "--window", str(window),
        "--min-id", str(min_id),
        "--min-ovlp", str(min_ovlp),
        "--top-out", str(top_out),
        "--expected-size", str(expected_size),
        "--size-tolerance", str(size_tolerance),
        "--max-depth", str(max_depth),
        "--max-cycles", str(max_cycles),
        "--max-linear", str(max_linear),
        "--max-cycles-any", str(max_cycles_any),
        "--fasta-top", str(fasta_top),
        "--out-prefix", out_prefix,
    ]
    log.info("disentangle_graph command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)

    report_txt = f"{out_prefix}.report.txt"
    with open(report_txt, "a") as fh:
        if proc.stderr:
            fh.write("\n[stderr]\n")
            fh.write(proc.stderr)

    if proc.returncode != 0:
        raise RuntimeError(
            f"disentangle_graph failed (exit={proc.returncode}). "
            f"See {report_txt} for details."
        )

    return {
        "gfa": gfa_path,
        "out_prefix": out_prefix,
        "terminal_overlaps_tsv": f"{out_prefix}.terminal_overlaps.tsv",
        "paths_tsv": f"{out_prefix}.paths.tsv",
        "cycles_any_tsv": f"{out_prefix}.cycles_any.tsv",
        "self_terminal_tsv": f"{out_prefix}.self_terminal.tsv",
        "path_sequences_fasta": f"{out_prefix}.path_sequences.fasta",
        "report_txt": report_txt,
    }


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------

def run_pipeline(
    reads: str,
    reference: str,
    outdir: str,
    threads: int = 4,
    # PAF filter options
    paf_len_a: int = 5000,
    paf_id_a: float = 0.90,
    paf_len_b: int = 3000,
    paf_id_b: float = 0.98,
    # Downsampler options
    target_coverage: int = 100,
    bin_size: int = 500,
    spanning_fraction: float = 0.80,
    # Assembler options
    run_assembler: bool = False,
    assembler: str = "flye",
    genome_size: str = "500k",
    # Disentangle options (always run when assembly is enabled)
    disentangle_window: int = 500,
    disentangle_min_id: float = 75.0,
    disentangle_min_ovlp: int = 100,
    disentangle_top_out: int = 5,
    disentangle_expected_size: int = 494000,
    disentangle_size_tolerance: float = 0.20,
    disentangle_max_depth: int = 8,
    disentangle_max_cycles: int = 30,
    disentangle_max_linear: int = 50,
    disentangle_max_cycles_any: int = 200,
    disentangle_fasta_top: int = 20,
    # Reporting
    report_name: str = "report.html",
) -> None:
    start_time = datetime.now()

    # Resolve chromosome names from reference
    chrom_names = _get_chrom_names(reference)
    if not chrom_names:
        raise ValueError(f"No sequences found in reference: {reference}")
    log.info("Reference chromosomes (%d): %s", len(chrom_names), chrom_names)

    # Create output directories
    dir_per_chrom = Path(outdir) / "01_per_chrom"
    dir_asm = Path(outdir) / "02_assembly"
    dir_report = Path(outdir) / "report"

    for d in [dir_per_chrom, dir_report]:
        d.mkdir(parents=True, exist_ok=True)

    report_path = str(dir_report / report_name)
    params_path = str(Path(outdir) / "pipeline_params.json")

    # Record parameters
    params = {
        "reads": os.path.abspath(reads),
        "reference": os.path.abspath(reference),
        "outdir": os.path.abspath(outdir),
        "threads": threads,
        "chromosomes": chrom_names,
        "paf_len_a": paf_len_a,
        "paf_id_a": paf_id_a,
        "paf_len_b": paf_len_b,
        "paf_id_b": paf_id_b,
        "target_coverage": target_coverage,
        "bin_size": bin_size,
        "spanning_fraction": spanning_fraction,
        "run_assembler": run_assembler,
        "assembler": assembler if run_assembler else "N/A",
        "genome_size": genome_size if run_assembler else "N/A",
        "disentangle_expected_size": disentangle_expected_size,
        "disentangle_size_tolerance": disentangle_size_tolerance,
        "started_at": start_time.isoformat(),
    }
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)
    log.info("Parameters saved to %s", params_path)

    # -----------------------------------------------------------------------
    # Step 1: Per-chromosome PAF filter
    # -----------------------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Per-chromosome PAF filter")
    log.info("=" * 60)

    chrom_filtered = _run_paf_filter_per_chrom(
        reads=reads,
        reference=reference,
        chrom_names=chrom_names,
        out_dir=dir_per_chrom,
        threads=threads,
        len_a=paf_len_a,
        id_a=paf_id_a,
        len_b=paf_len_b,
        id_b=paf_id_b,
    )

    # -----------------------------------------------------------------------
    # Step 2: Per-chromosome downsampler
    # -----------------------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 2 — Per-chromosome downsampler")
    log.info("=" * 60)

    chrom_downsampled: Dict[str, str] = {}
    all_ds_stats: Dict[str, Dict] = {}

    for chrom in chrom_names:
        filtered_fq = chrom_filtered[chrom]
        chrom_dir = dir_per_chrom / chrom
        ds_fastq = str(chrom_dir / "downsampled.fastq")
        ds_stats_path = str(chrom_dir / "downsample_stats.json")

        stats = _downsample_chrom(
            filtered_reads=filtered_fq,
            reference=reference,
            chrom_name=chrom,
            out_fastq=ds_fastq,
            stats_path=ds_stats_path,
            threads=threads,
            target_coverage=target_coverage,
            bin_size=bin_size,
            spanning_fraction=spanning_fraction,
        )
        chrom_downsampled[chrom] = ds_fastq
        all_ds_stats[chrom] = stats

    # -----------------------------------------------------------------------
    # Step 3 (optional): Per-chromosome assembly + mandatory disentangle
    # -----------------------------------------------------------------------
    chrom_assembly_paths: Dict[str, Optional[str]] = {}
    chrom_disentangle_results: Dict[str, Optional[Dict]] = {}

    if run_assembler:
        fn = ASSEMBLERS.get(assembler)
        if fn is None:
            log.warning("Unknown assembler '%s'. Skipping.", assembler)
        else:
            for chrom in chrom_names:
                log.info("=" * 60)
                log.info("STEP 3 — Assembly of %s with %s", chrom, assembler)
                log.info("=" * 60)

                chrom_asm_dir = dir_asm / chrom
                chrom_asm_dir.mkdir(parents=True, exist_ok=True)

                ds_reads = chrom_downsampled[chrom]
                assembly_path: Optional[str] = None

                try:
                    assembly_path = fn(
                        reads=ds_reads,
                        outdir=str(chrom_asm_dir),
                        threads=threads,
                        genome_size=genome_size,
                    )
                    log.info("Assembly for %s written to %s", chrom, assembly_path)
                except Exception as exc:
                    log.error("Assembly of %s failed: %s", chrom, exc)

                chrom_assembly_paths[chrom] = assembly_path

                # Disentangle is MANDATORY after every successful assembly
                gfa_path = str(chrom_asm_dir / "assembly_graph.gfa")
                if not os.path.exists(gfa_path):
                    log.warning(
                        "disentangle skipped for %s: GFA not found at %s",
                        chrom, gfa_path,
                    )
                    chrom_disentangle_results[chrom] = None
                    continue

                log.info("=" * 60)
                log.info("STEP 3b — Disentangle graph for %s (mandatory)", chrom)
                log.info("=" * 60)
                dis_prefix = str(chrom_asm_dir / "disentangle")
                try:
                    dis_result = _run_disentangle(
                        gfa_path=gfa_path,
                        out_prefix=dis_prefix,
                        window=disentangle_window,
                        min_id=disentangle_min_id,
                        min_ovlp=disentangle_min_ovlp,
                        top_out=disentangle_top_out,
                        expected_size=disentangle_expected_size,
                        size_tolerance=disentangle_size_tolerance,
                        max_depth=disentangle_max_depth,
                        max_cycles=disentangle_max_cycles,
                        max_linear=disentangle_max_linear,
                        max_cycles_any=disentangle_max_cycles_any,
                        fasta_top=disentangle_fasta_top,
                    )
                    log.info("disentangle_graph for %s: %s", chrom, dis_result["report_txt"])
                    chrom_disentangle_results[chrom] = dis_result
                except Exception as exc:
                    log.error("disentangle_graph for %s failed: %s", chrom, exc)
                    chrom_disentangle_results[chrom] = None
    else:
        log.info("Assembly step skipped (use --run-assembler to enable)")

    # -----------------------------------------------------------------------
    # Step 4: HTML Report
    # -----------------------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 4 — Generating HTML Report")
    log.info("=" * 60)

    params["finished_at"] = datetime.now().isoformat()
    elapsed = (datetime.now() - start_time).total_seconds()
    params["elapsed_seconds"] = round(elapsed, 1)
    if chrom_disentangle_results:
        params["disentangle_results"] = {
            k: v for k, v in chrom_disentangle_results.items() if v is not None
        }
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)

    # Use the first non-empty downsampled file as "final_reads" for report
    final_reads_for_report = next(
        (p for p in chrom_downsampled.values() if os.path.exists(p) and os.path.getsize(p) > 0),
        reads,
    )
    # Aggregate ds_stats across chromosomes for report
    total_output_reads = sum(s["output_reads"] for s in all_ds_stats.values())
    agg_ds_stats: Dict = {
        "input_alignments": sum(s["input_alignments"] for s in all_ds_stats.values()),
        "spanning_reads": sum(s["spanning_reads"] for s in all_ds_stats.values()),
        "selected_reads": sum(s["selected_reads"] for s in all_ds_stats.values()),
        "output_reads": total_output_reads,
        "target_coverage": target_coverage,
        "achieved_coverage": {
            chrom: all_ds_stats[chrom]["achieved_coverage"]
            for chrom in chrom_names
        },
        "per_chrom": all_ds_stats,
    }

    report_mod.generate_report(
        raw_reads=reads,
        scrubbed_reads=reads,       # no scrub step; raw == scrubbed
        final_reads=final_reads_for_report,
        output_html=report_path,
        numt_stats={"filter_mode": "paf_per_chrom", "scrub_removed": False},
        downsample_stats=agg_ds_stats,
        assembly_path=next(
            (p for p in chrom_assembly_paths.values() if p is not None), None
        ),
        parameters={
            "Threads": threads,
            "Chromosomes": ", ".join(chrom_names),
            "PAF len_A / id_A": f"{paf_len_a} / {paf_id_a}",
            "PAF len_B / id_B": f"{paf_len_b} / {paf_id_b}",
            "Target coverage": f"{target_coverage}x",
            "Bin size (bp)": bin_size,
            "Spanning fraction": spanning_fraction,
            "Assembler": assembler if run_assembler else "Not run",
            "Genome size estimate": genome_size if run_assembler else "N/A",
            "Elapsed time (s)": params["elapsed_seconds"],
        },
    )

    log.info("=" * 60)
    log.info("Pipeline complete in %.1f seconds", elapsed)
    log.info("Report: %s", report_path)
    log.info("=" * 60)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "MitoAssembler: End-to-end plant mitogenome assembly pipeline.\n"
            "Per-chromosome PAF filter → per-chromosome downsampling → "
            "(optional) per-chromosome assembly → disentangle graph → HTML report."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    p.add_argument("--reads", required=True, help="Raw ONT reads (FASTQ/FASTA)")
    p.add_argument("--reference", required=True, help="Mitochondrial reference FASTA")
    p.add_argument("--outdir", required=True, help="Output directory (created if absent)")

    # General
    p.add_argument("--threads", type=int, default=4)

    # PAF filter
    grp_paf = p.add_argument_group("Per-chromosome PAF filter")
    grp_paf.add_argument("--paf-len-a", type=int, default=5000, dest="paf_len_a")
    grp_paf.add_argument("--paf-id-a", type=float, default=0.90, dest="paf_id_a")
    grp_paf.add_argument("--paf-len-b", type=int, default=3000, dest="paf_len_b")
    grp_paf.add_argument("--paf-id-b", type=float, default=0.98, dest="paf_id_b")

    # Downsampler
    grp_ds = p.add_argument_group("Per-chromosome Downsampler")
    grp_ds.add_argument("--target-coverage", type=int, default=100, dest="target_coverage")
    grp_ds.add_argument("--bin-size", type=int, default=500, dest="bin_size")
    grp_ds.add_argument("--spanning-fraction", type=float, default=0.80,
                        dest="spanning_fraction")

    # Assembler
    grp_asm = p.add_argument_group("Assembler (optional — per chromosome)")
    grp_asm.add_argument("--run-assembler", action="store_true", dest="run_assembler",
                         help="Run per-chromosome assembly after downsampling")
    grp_asm.add_argument("--assembler", default="flye",
                         choices=list(ASSEMBLERS.keys()))
    grp_asm.add_argument("--genome-size", default="500k", dest="genome_size",
                         help="Estimated genome size for assembler (e.g. 500k)")

    # Disentangle (mandatory when assembly runs)
    grp_dis = p.add_argument_group("Disentangle graph (runs automatically with --run-assembler)")
    grp_dis.add_argument("--disentangle-window", type=int, default=500, dest="disentangle_window")
    grp_dis.add_argument("--disentangle-min-id", type=float, default=75.0, dest="disentangle_min_id")
    grp_dis.add_argument("--disentangle-min-ovlp", type=int, default=100, dest="disentangle_min_ovlp")
    grp_dis.add_argument("--disentangle-top-out", type=int, default=5, dest="disentangle_top_out")
    grp_dis.add_argument("--disentangle-expected-size", type=int, default=494000,
                         dest="disentangle_expected_size")
    grp_dis.add_argument("--disentangle-size-tolerance", type=float, default=0.20,
                         dest="disentangle_size_tolerance")
    grp_dis.add_argument("--disentangle-max-depth", type=int, default=8, dest="disentangle_max_depth")
    grp_dis.add_argument("--disentangle-max-cycles", type=int, default=30, dest="disentangle_max_cycles")
    grp_dis.add_argument("--disentangle-max-linear", type=int, default=50, dest="disentangle_max_linear")
    grp_dis.add_argument("--disentangle-max-cycles-any", type=int, default=200,
                         dest="disentangle_max_cycles_any")
    grp_dis.add_argument("--disentangle-fasta-top", type=int, default=20, dest="disentangle_fasta_top")

    # Report
    p.add_argument("--report-name", default="report.html", dest="report_name")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_pipeline(
        reads=args.reads,
        reference=args.reference,
        outdir=args.outdir,
        threads=args.threads,
        paf_len_a=args.paf_len_a,
        paf_id_a=args.paf_id_a,
        paf_len_b=args.paf_len_b,
        paf_id_b=args.paf_id_b,
        target_coverage=args.target_coverage,
        bin_size=args.bin_size,
        spanning_fraction=args.spanning_fraction,
        run_assembler=args.run_assembler,
        assembler=args.assembler,
        genome_size=args.genome_size,
        disentangle_window=args.disentangle_window,
        disentangle_min_id=args.disentangle_min_id,
        disentangle_min_ovlp=args.disentangle_min_ovlp,
        disentangle_top_out=args.disentangle_top_out,
        disentangle_expected_size=args.disentangle_expected_size,
        disentangle_size_tolerance=args.disentangle_size_tolerance,
        disentangle_max_depth=args.disentangle_max_depth,
        disentangle_max_cycles=args.disentangle_max_cycles,
        disentangle_max_linear=args.disentangle_max_linear,
        disentangle_max_cycles_any=args.disentangle_max_cycles_any,
        disentangle_fasta_top=args.disentangle_fasta_top,
        report_name=args.report_name,
    )


if __name__ == "__main__":
    main()
