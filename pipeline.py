#!/usr/bin/env python3
"""
pipeline.py — MitoAssembler End-to-End Pipeline
================================================
Single-command wrapper that orchestrates:
  1. NUMT Scrubber    (numt_scrubber.py)
  2. Uniform Downsampler (downsampler.py)
  3. HTML Report      (report.py)

Optionally runs an assembler (default: Flye in --nano-hq mode) on the
final clean reads to produce a mitogenome assembly.

Usage
-----
    python pipeline.py \\
        --reads raw_reads.fastq \\
        --reference mito_ref.fasta \\
        --outdir results/ \\
        [--threads 16] \\
        [--target-coverage 100] \\
        [--run-assembler] \\
        [--assembler flye]  \\
        [--genome-size 500k]

Output layout (inside --outdir)
--------------------------------
    01_numt_scrubbed/
        scrubbed.fastq
        numt_stats.json
    02_downsampled/
        downsampled.fastq
        downsample_stats.json
    03_assembly/           (only when --run-assembler)
        assembly.fasta
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
from typing import Dict, Optional, Set

import pysam

# Local modules (same directory)
sys.path.insert(0, os.path.dirname(__file__))

import numt_scrubber
import downsampler
import report as report_mod

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# PAF dual-threshold filter  (replicates old_methods/align2.sh)
# ---------------------------------------------------------------------------

def _run_paf_filter(
    reads: str,
    reference: str,
    output_fastq: str,
    threads: int,
    len_a: int = 5000,
    id_a: float = 0.90,
    len_b: int = 3000,
    id_b: float = 0.98,
) -> Dict:
    """
    Dual-threshold PAF filter that replicates old_methods/align2.sh.

    Keep a read when EITHER condition is true:
      A) aligned_length >= len_a  AND  identity >= id_a  (long reads, relaxed id)
      B) aligned_length >= len_b  AND  identity >= id_b  (medium reads, strict id)

    identity = matching_bases / aligned_length  (PAF col10 / col11)
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

        keep_names: Set[str] = set()
        total_aln = 0
        passed_a = 0
        passed_b = 0

        with open(paf_path) as paf_fh:
            for line in paf_fh:
                cols = line.split("\t")
                if len(cols) < 12:
                    continue
                total_aln += 1
                read_name = cols[0]
                aln_len = int(cols[10])
                matches = int(cols[9])
                if aln_len == 0:
                    continue
                identity = matches / aln_len
                if aln_len >= len_a and identity >= id_a:
                    keep_names.add(read_name)
                    passed_a += 1
                elif aln_len >= len_b and identity >= id_b:
                    keep_names.add(read_name)
                    passed_b += 1

        log.info(
            "PAF filter: %d alignments → %d unique reads (cond_A=%d cond_B=%d)",
            total_aln, len(keep_names), passed_a, passed_b,
        )

    written = 0
    with pysam.FastxFile(reads) as fin, open(output_fastq, "w") as fout:
        for entry in fin:
            if entry.name in keep_names:
                qual = entry.quality or ""
                fout.write(f"@{entry.name}\n{entry.sequence}\n+\n{qual}\n")
                written += 1

    log.info("PAF filter: wrote %d reads to %s", written, output_fastq)

    return {
        "filter_mode": "paf",
        "total_alignments": total_aln,
        "passed_condition_a": passed_a,
        "passed_condition_b": passed_b,
        "unique_reads_kept": len(keep_names),
        "output_reads": written,
        "thresholds": {"len_a": len_a, "id_a": id_a, "len_b": len_b, "id_b": id_b},
    }

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


def _run_check_joins(
    gfa_path: str,
    out_tsv: str,
    summary_txt: str,
    window: int,
    min_id: float,
    min_ovlp: int,
    expected_size: int,
    size_tolerance: float,
) -> Dict:
    check_joins_script = os.path.join(os.path.dirname(__file__), "check_joins.py")
    if not os.path.exists(check_joins_script):
        raise FileNotFoundError(f"check_joins.py not found at {check_joins_script}")

    cmd = [
        sys.executable,
        check_joins_script,
        gfa_path,
        "--window",
        str(window),
        "--min-id",
        str(min_id),
        "--min-ovlp",
        str(min_ovlp),
        "--expected-size",
        str(expected_size),
        "--size-tolerance",
        str(size_tolerance),
        "--out",
        out_tsv,
    ]
    log.info("check_joins command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)

    with open(summary_txt, "w") as fh:
        if proc.stdout:
            fh.write(proc.stdout)
        if proc.stderr:
            fh.write("\n[stderr]\n")
            fh.write(proc.stderr)

    if proc.returncode != 0:
        raise RuntimeError(
            f"check_joins failed (exit={proc.returncode}). "
            f"See {summary_txt} for details."
        )

    return {
        "gfa": gfa_path,
        "candidate_joins_tsv": out_tsv,
        "summary_txt": summary_txt,
    }


ASSEMBLERS = {
    "flye": _run_flye,
}


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------

def run_pipeline(
    reads: str,
    reference: str,
    outdir: str,
    threads: int = 4,
    # Filter mode
    filter_mode: str = "numt",
    # NUMT scrubber options (filter_mode='numt')
    min_mapq: int = 20,
    min_mapped_fraction: float = 0.80,
    max_softclip_fraction: float = 0.20,
    coverage_window: int = 500,
    coverage_threshold_ratio: float = 0.20,
    # PAF dual-threshold options (filter_mode='paf')
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
    run_check_joins: bool = False,
    joins_window: int = 500,
    joins_min_id: float = 75.0,
    joins_min_ovlp: int = 100,
    joins_expected_size: int = 494000,
    joins_size_tolerance: float = 0.20,
    # Reporting
    report_name: str = "report.html",
) -> None:
    start_time = datetime.now()

    # Create output directories
    dir_scrubbed = Path(outdir) / "01_numt_scrubbed"
    dir_down = Path(outdir) / "02_downsampled"
    dir_asm = Path(outdir) / "03_assembly"
    dir_report = Path(outdir) / "report"

    for d in [dir_scrubbed, dir_down, dir_report]:
        d.mkdir(parents=True, exist_ok=True)

    scrubbed_reads = str(dir_scrubbed / "scrubbed.fastq")
    numt_stats_path = str(dir_scrubbed / "numt_stats.json")
    downsampled_reads = str(dir_down / "downsampled.fastq")
    ds_stats_path = str(dir_down / "downsample_stats.json")
    report_path = str(dir_report / report_name)
    params_path = str(Path(outdir) / "pipeline_params.json")

    # Record parameters
    params = {
        "reads": os.path.abspath(reads),
        "reference": os.path.abspath(reference),
        "outdir": os.path.abspath(outdir),
        "threads": threads,
        "min_mapq": min_mapq,
        "min_mapped_fraction": min_mapped_fraction,
        "max_softclip_fraction": max_softclip_fraction,
        "coverage_window": coverage_window,
        "coverage_threshold_ratio": coverage_threshold_ratio,
        "target_coverage": target_coverage,
        "bin_size": bin_size,
        "spanning_fraction": spanning_fraction,
        "run_assembler": run_assembler,
        "assembler": assembler if run_assembler else "N/A",
        "genome_size": genome_size if run_assembler else "N/A",
        "started_at": start_time.isoformat(),
        "filter_mode": filter_mode,
        "paf_len_a": paf_len_a,
        "paf_id_a": paf_id_a,
        "paf_len_b": paf_len_b,
        "paf_id_b": paf_id_b,
        "run_check_joins": run_check_joins,
        "joins_window": joins_window,
        "joins_min_id": joins_min_id,
        "joins_min_ovlp": joins_min_ovlp,
        "joins_expected_size": joins_expected_size,
        "joins_size_tolerance": joins_size_tolerance,
    }
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)
    log.info("Parameters saved to %s", params_path)

    # -----------------------------------------------------------------------
    # Step 1: Read Filtering
    # -----------------------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Read Filtering (mode: %s)", filter_mode)
    log.info("=" * 60)

    if filter_mode == "paf":
        # PAF dual-threshold filter — replicates old_methods/align2.sh
        # Bypasses BAM-based NUMT scrubber; uses raw PAF identity thresholds.
        numt_stats = _run_paf_filter(
            reads=reads,
            reference=reference,
            output_fastq=scrubbed_reads,
            threads=threads,
            len_a=paf_len_a,
            id_a=paf_id_a,
            len_b=paf_len_b,
            id_b=paf_id_b,
        )
        written = numt_stats["output_reads"]
    else:
        # Default: BAM-based NUMT scrubber with per-contig coverage baseline
        with tempfile.TemporaryDirectory() as tmpdir:
            bam_path = os.path.join(tmpdir, "numt_aligned.bam")

            numt_scrubber.run_minimap2(
                reference=reference,
                reads=reads,
                bam_out=bam_path,
                threads=threads,
            )

            cov_data = numt_scrubber.compute_sliding_window_coverage(
                bam_path, window=coverage_window
            )
            per_contig_baselines = numt_scrubber.estimate_baselines(cov_data)
            for contig, bl in per_contig_baselines.items():
                log.info("  Contig %-30s  baseline=%.1fx", contig, bl)

            clean_names, numt_stats = numt_scrubber.classify_reads(
                bam_path=bam_path,
                per_contig_baselines=per_contig_baselines,
                coverage_threshold_ratio=coverage_threshold_ratio,
                min_mapq=min_mapq,
                min_mapped_fraction=min_mapped_fraction,
                max_softclip_fraction=max_softclip_fraction,
            )

            written = numt_scrubber.extract_clean_reads(reads, clean_names, scrubbed_reads)
            log.info("NUMT scrubber: %d → %d reads", numt_stats["total_alignments"], written)
            numt_stats["output_reads"] = written

    with open(numt_stats_path, "w") as fh:
        json.dump(numt_stats, fh, indent=2)

    # -----------------------------------------------------------------------
    # Step 2: Downsampler
    # -----------------------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 2 — Uniform Downsampler")
    log.info("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = os.path.join(tmpdir, "ds_aligned.bam")

        downsampler.run_minimap2(
            reference=reference,
            reads=scrubbed_reads,
            bam_out=bam_path,
            threads=threads,
        )

        ds_reads = downsampler.load_read_info(bam_path, min_mapq=10)
        bam = pysam.AlignmentFile(bam_path, "rb")
        ref_lengths = dict(zip(bam.references, bam.lengths))
        bam.close()

        spanning = downsampler.find_spanning_reads(ds_reads, ref_lengths, spanning_fraction)
        selected = downsampler.greedy_bin_selection(
            ds_reads, ref_lengths, target_coverage=target_coverage, bin_size=bin_size
        )
        keep = selected | spanning

        ds_written = downsampler.extract_reads(scrubbed_reads, keep, downsampled_reads)
        log.info("Downsampler: %d → %d reads", len(ds_reads), ds_written)

        achieved = downsampler.compute_achieved_coverage(keep, ds_reads, ref_lengths)

    ds_stats: Dict = {
        "input_alignments": len(ds_reads),
        "spanning_reads": len(spanning),
        "selected_reads": len(keep),
        "output_reads": ds_written,
        "target_coverage": target_coverage,
        "achieved_coverage": achieved,
        "ref_lengths": ref_lengths,  # needed by report for length-weighted coverage mean
    }
    with open(ds_stats_path, "w") as fh:
        json.dump(ds_stats, fh, indent=2)

    # -----------------------------------------------------------------------
    # Step 3 (optional): Assembly
    # -----------------------------------------------------------------------
    assembly_path: Optional[str] = None
    join_results: Optional[Dict] = None
    if run_assembler:
        log.info("=" * 60)
        log.info("STEP 3 — Assembly with %s", assembler)
        log.info("=" * 60)
        dir_asm.mkdir(parents=True, exist_ok=True)
        fn = ASSEMBLERS.get(assembler)
        if fn is None:
            log.warning("Unknown assembler '%s'. Skipping.", assembler)
        else:
            try:
                assembly_path = fn(
                    reads=downsampled_reads,
                    outdir=str(dir_asm),
                    threads=threads,
                    genome_size=genome_size,
                )
                log.info("Assembly written to %s", assembly_path)
            except Exception as exc:
                log.error("Assembly failed: %s", exc)

        if run_check_joins:
            gfa_path = os.path.join(str(dir_asm), "assembly_graph.gfa")
            if not os.path.exists(gfa_path):
                log.warning("check_joins skipped: GFA not found at %s", gfa_path)
            else:
                log.info("=" * 60)
                log.info("STEP 3b — Graph Join Analysis (check_joins)")
                log.info("=" * 60)
                joins_tsv = os.path.join(str(dir_asm), "candidate_joins.tsv")
                joins_summary = os.path.join(str(dir_asm), "check_joins_summary.txt")
                try:
                    join_results = _run_check_joins(
                        gfa_path=gfa_path,
                        out_tsv=joins_tsv,
                        summary_txt=joins_summary,
                        window=joins_window,
                        min_id=joins_min_id,
                        min_ovlp=joins_min_ovlp,
                        expected_size=joins_expected_size,
                        size_tolerance=joins_size_tolerance,
                    )
                    log.info("check_joins results: %s", join_results["candidate_joins_tsv"])
                except Exception as exc:
                    log.error("check_joins failed: %s", exc)
    else:
        log.info("Assembly step skipped (use --run-assembler to enable)")

    if join_results is not None:
        params["check_joins_results"] = join_results

    # -----------------------------------------------------------------------
    # Step 4: HTML Report
    # -----------------------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 4 — Generating HTML Report")
    log.info("=" * 60)

    params["finished_at"] = datetime.now().isoformat()
    elapsed = (datetime.now() - start_time).total_seconds()
    params["elapsed_seconds"] = round(elapsed, 1)
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)

    report_mod.generate_report(
        raw_reads=reads,
        scrubbed_reads=scrubbed_reads,
        final_reads=downsampled_reads,
        output_html=report_path,
        numt_stats=numt_stats,
        downsample_stats=ds_stats,
        assembly_path=assembly_path,
        parameters={
            "Threads": threads,
            "Min MAPQ": min_mapq,
            "Min mapped fraction": min_mapped_fraction,
            "Max softclip fraction": max_softclip_fraction,
            "Coverage window (bp)": coverage_window,
            "Coverage threshold ratio": coverage_threshold_ratio,
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
            "Runs NUMT scrubbing → uniform downsampling → (optional) assembly → HTML report."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    p.add_argument("--reads", required=True, help="Raw ONT reads (FASTQ/FASTA)")
    p.add_argument("--reference", required=True, help="Mitochondrial reference FASTA")
    p.add_argument("--outdir", required=True, help="Output directory (created if absent)")

    # General
    p.add_argument("--threads", type=int, default=4)

    # Filter mode
    grp_filter = p.add_argument_group("Read Filtering Mode")
    grp_filter.add_argument(
        "--filter-mode",
        choices=["numt", "paf"],
        default="numt",
        dest="filter_mode",
        help=(
            "Filtering strategy: 'numt' uses BAM+coverage NUMT scrubber; "
            "'paf' uses old_methods dual-threshold PAF filtering"
        ),
    )
    grp_filter.add_argument("--paf-len-a", type=int, default=5000, dest="paf_len_a")
    grp_filter.add_argument("--paf-id-a", type=float, default=0.90, dest="paf_id_a")
    grp_filter.add_argument("--paf-len-b", type=int, default=3000, dest="paf_len_b")
    grp_filter.add_argument("--paf-id-b", type=float, default=0.98, dest="paf_id_b")

    # NUMT scrubber
    grp_numt = p.add_argument_group("NUMT Scrubber")
    grp_numt.add_argument("--min-mapq", type=int, default=20, dest="min_mapq")
    grp_numt.add_argument("--min-mapped-fraction", type=float, default=0.80,
                          dest="min_mapped_fraction")
    grp_numt.add_argument("--max-softclip-fraction", type=float, default=0.20,
                          dest="max_softclip_fraction")
    grp_numt.add_argument("--coverage-window", type=int, default=500,
                          dest="coverage_window")
    grp_numt.add_argument("--coverage-threshold-ratio", type=float, default=0.20,
                          dest="coverage_threshold_ratio")

    # Downsampler
    grp_ds = p.add_argument_group("Downsampler")
    grp_ds.add_argument("--target-coverage", type=int, default=100, dest="target_coverage")
    grp_ds.add_argument("--bin-size", type=int, default=500, dest="bin_size")
    grp_ds.add_argument("--spanning-fraction", type=float, default=0.80,
                        dest="spanning_fraction")

    # Assembler
    grp_asm = p.add_argument_group("Assembler (optional)")
    grp_asm.add_argument("--run-assembler", action="store_true", dest="run_assembler",
                         help="Run assembly after downsampling")
    grp_asm.add_argument("--assembler", default="flye",
                         choices=list(ASSEMBLERS.keys()))
    grp_asm.add_argument("--genome-size", default="500k", dest="genome_size",
                         help="Estimated genome size for assembler (e.g. 500k)")
    grp_asm.add_argument(
        "--run-check-joins",
        action="store_true",
        dest="run_check_joins",
        help="Run check_joins.py on Flye assembly_graph.gfa when assembly is enabled",
    )
    grp_asm.add_argument("--joins-window", type=int, default=500, dest="joins_window")
    grp_asm.add_argument("--joins-min-id", type=float, default=75.0, dest="joins_min_id")
    grp_asm.add_argument("--joins-min-ovlp", type=int, default=100, dest="joins_min_ovlp")
    grp_asm.add_argument(
        "--joins-expected-size",
        type=int,
        default=494000,
        dest="joins_expected_size",
    )
    grp_asm.add_argument(
        "--joins-size-tolerance",
        type=float,
        default=0.20,
        dest="joins_size_tolerance",
    )

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
        filter_mode=args.filter_mode,
        min_mapq=args.min_mapq,
        min_mapped_fraction=args.min_mapped_fraction,
        max_softclip_fraction=args.max_softclip_fraction,
        coverage_window=args.coverage_window,
        coverage_threshold_ratio=args.coverage_threshold_ratio,
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
        run_check_joins=args.run_check_joins,
        joins_window=args.joins_window,
        joins_min_id=args.joins_min_id,
        joins_min_ovlp=args.joins_min_ovlp,
        joins_expected_size=args.joins_expected_size,
        joins_size_tolerance=args.joins_size_tolerance,
        report_name=args.report_name,
    )


if __name__ == "__main__":
    main()
