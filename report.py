#!/usr/bin/env python3
"""
report.py — HTML QC Report Generator
======================================
Produces a self-contained HTML report summarising the MitoAssembler pipeline run.

Contents
--------
- Run metadata (date, parameters)
- Read counts at each pipeline stage (input → NUMT-scrubbed → downsampled)
- Per-stage Q-score distribution (histogram)
- Coverage plots before/after downsampling
- NUMT / chimeric read breakdown (stacked bar)
- Final assembly statistics (if a FASTA assembly is provided)

All charts use Chart.js (loaded from CDN) embedded in a single HTML file.
No external Python plotting dependencies are required.

Usage (standalone)
------------------
    python report.py \\
        --numt-stats numt_stats.json \\
        --downsample-stats downsample_stats.json \\
        --raw-reads raw.fastq \\
        --scrubbed-reads scrubbed.fastq \\
        --final-reads final.fastq \\
        [--assembly assembly.fasta] \\
        --output report.html

It is also importable — call generate_report(**kwargs) directly from pipeline.py.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import pysam


# ---------------------------------------------------------------------------
# Read statistics helpers
# ---------------------------------------------------------------------------

def read_stats(path: str) -> Dict[str, Any]:
    """Compute basic statistics from a FASTQ/FASTA file using pysam."""
    lengths: List[int] = []
    quals: List[float] = []
    with pysam.FastxFile(path) as fh:
        for entry in fh:
            lengths.append(len(entry.sequence))
            if entry.quality:
                q = [ord(c) - 33 for c in entry.quality]
                quals.append(sum(q) / len(q))

    if not lengths:
        return {
            "count": 0, "total_bases": 0, "mean_length": 0, "n50": 0,
            "mean_q": 0, "q_histogram": [],
        }

    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    cumsum = 0
    n50 = 0
    for l in lengths:
        cumsum += l
        if cumsum >= half:
            n50 = l
            break

    # Q-score histogram bins: 0–5, 5–10, …, 40–45
    q_bins = list(range(0, 46, 5))
    hist = [0] * (len(q_bins) - 1)
    for q in quals:
        idx = min(int(q // 5), len(hist) - 1)
        hist[idx] += 1

    return {
        "count": len(lengths),
        "total_bases": total,
        "mean_length": int(total / len(lengths)),
        "n50": n50,
        "mean_q": round(sum(quals) / len(quals), 2) if quals else 0,
        "q_histogram": hist,
        "q_bins": [f"Q{q_bins[i]}–Q{q_bins[i+1]}" for i in range(len(q_bins) - 1)],
    }


def fasta_assembly_stats(path: str) -> Dict[str, Any]:
    """Return contig-level statistics for a FASTA assembly."""
    contigs: List[Tuple[str, int]] = []
    with pysam.FastxFile(path) as fh:
        for entry in fh:
            contigs.append((entry.name, len(entry.sequence)))

    if not contigs:
        return {"n_contigs": 0, "total_length": 0, "longest": 0, "contigs": []}

    lengths = [l for _, l in contigs]
    return {
        "n_contigs": len(contigs),
        "total_length": sum(lengths),
        "longest": max(lengths),
        "contigs": [{"name": n, "length": l} for n, l in contigs],
    }


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>MitoAssembler Pipeline Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
  :root {{
    --bg: #f9fafb;
    --card: #ffffff;
    --accent: #2563eb;
    --text: #1f2937;
    --sub: #6b7280;
    --border: #e5e7eb;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: 'Segoe UI', system-ui, sans-serif;
    background: var(--bg);
    color: var(--text);
    padding: 2rem;
  }}
  h1 {{ font-size: 1.8rem; margin-bottom: 0.25rem; color: var(--accent); }}
  h2 {{ font-size: 1.2rem; margin: 1.5rem 0 0.75rem; border-bottom: 2px solid var(--accent); padding-bottom: 0.3rem; }}
  h3 {{ font-size: 1rem; margin-bottom: 0.5rem; color: var(--sub); }}
  .meta {{ font-size: 0.85rem; color: var(--sub); margin-bottom: 2rem; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 1rem; margin-bottom: 1.5rem; }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.2rem;
    box-shadow: 0 1px 3px rgba(0,0,0,.06);
  }}
  .card .value {{ font-size: 2rem; font-weight: 700; color: var(--accent); }}
  .card .label {{ font-size: 0.8rem; color: var(--sub); margin-top: 0.2rem; }}
  .chart-wrap {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,.06);
  }}
  canvas {{ max-height: 320px; }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.9rem;
    background: var(--card);
    border-radius: 10px;
    overflow: hidden;
    box-shadow: 0 1px 3px rgba(0,0,0,.06);
    margin-bottom: 1.5rem;
  }}
  th, td {{ padding: 0.7rem 1rem; border-bottom: 1px solid var(--border); text-align: left; }}
  th {{ background: var(--accent); color: white; }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: #f0f4ff; }}
  .badge {{
    display: inline-block;
    padding: 0.2rem 0.6rem;
    border-radius: 999px;
    font-size: 0.75rem;
    font-weight: 600;
  }}
  .badge-green {{ background: #d1fae5; color: #065f46; }}
  .badge-yellow {{ background: #fef3c7; color: #92400e; }}
  .badge-red {{ background: #fee2e2; color: #991b1b; }}
  footer {{ text-align: center; font-size: 0.8rem; color: var(--sub); margin-top: 3rem; }}
</style>
</head>
<body>

<h1>MitoAssembler Pipeline Report</h1>
<p class="meta">Generated: {generated_at} &nbsp;|&nbsp; Pipeline: NUMT Scrubber &rarr; Uniform Downsampler</p>

<h2>Read Counts Across Pipeline Stages</h2>
<div class="grid">
  <div class="card"><div class="value">{raw_count}</div><div class="label">Raw input reads</div></div>
  <div class="card"><div class="value">{scrubbed_count}</div><div class="label">After NUMT scrubbing</div></div>
  <div class="card"><div class="value">{final_count}</div><div class="label">After downsampling</div></div>
  <div class="card"><div class="value">{numt_removed}</div><div class="label">NUMT / chimeric removed</div></div>
  <div class="card"><div class="value">{reduction_pct}%</div><div class="label">Total read reduction</div></div>
  <div class="card"><div class="value">{achieved_cov}</div><div class="label">Achieved mean coverage</div></div>
</div>

<div class="chart-wrap">
  <h3>Read Count per Stage</h3>
  <canvas id="stageBar"></canvas>
</div>

<h2>NUMT Scrubber Breakdown</h2>
<div class="chart-wrap">
  <h3>Reads removed by filter type</h3>
  <canvas id="numtPie"></canvas>
</div>

<h2>Q-Score Distribution</h2>
<div class="chart-wrap">
  <h3>Mean Q-score histogram (reads)</h3>
  <canvas id="qHist"></canvas>
</div>

<h2>Read Length Statistics</h2>
<table>
  <thead>
    <tr><th>Stage</th><th>Read Count</th><th>Total Bases (Mb)</th><th>Mean Length (bp)</th><th>N50 (bp)</th><th>Mean Q</th></tr>
  </thead>
  <tbody>
    <tr>
      <td>Raw input</td>
      <td>{raw_count}</td>
      <td>{raw_total_mb}</td>
      <td>{raw_mean_len:,}</td>
      <td>{raw_n50:,}</td>
      <td>{raw_mean_q}</td>
    </tr>
    <tr>
      <td>NUMT scrubbed</td>
      <td>{scrubbed_count}</td>
      <td>{scrubbed_total_mb}</td>
      <td>{scrubbed_mean_len:,}</td>
      <td>{scrubbed_n50:,}</td>
      <td>{scrubbed_mean_q}</td>
    </tr>
    <tr>
      <td>Downsampled (final)</td>
      <td>{final_count}</td>
      <td>{final_total_mb}</td>
      <td>{final_mean_len:,}</td>
      <td>{final_n50:,}</td>
      <td>{final_mean_q}</td>
    </tr>
  </tbody>
</table>

{assembly_section}

<h2>Per-Chromosome Coverage After Downsampling</h2>
<table>
  <thead><tr><th>Chromosome / Contig</th><th>Achieved Coverage</th></tr></thead>
  <tbody>
    {per_contig_cov_rows}
  </tbody>
</table>

<h2>Pipeline Parameters</h2>
<table>
  <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
  <tbody>
    {param_rows}
  </tbody>
</table>

<footer>MitoAssembler &mdash; Plant Mitogenome Assembly Pipeline</footer>

<script>
// Stage bar chart
new Chart(document.getElementById('stageBar'), {{
  type: 'bar',
  data: {{
    labels: ['Raw input', 'NUMT scrubbed', 'Downsampled'],
    datasets: [{{
      label: 'Read count',
      data: [{raw_count}, {scrubbed_count}, {final_count}],
      backgroundColor: ['#93c5fd', '#34d399', '#f59e0b'],
      borderRadius: 6,
    }}]
  }},
  options: {{
    plugins: {{ legend: {{ display: false }} }},
    scales: {{ y: {{ beginAtZero: true }} }}
  }}
}});

// NUMT pie chart
new Chart(document.getElementById('numtPie'), {{
  type: 'doughnut',
  data: {{
    labels: {numt_labels_json},
    datasets: [{{
      data: {numt_values_json},
      backgroundColor: ['#34d399', '#f87171', '#fb923c', '#a78bfa', '#38bdf8'],
    }}]
  }},
  options: {{
    plugins: {{
      legend: {{ position: 'right' }}
    }}
  }}
}});

// Q-score histogram
new Chart(document.getElementById('qHist'), {{
  type: 'bar',
  data: {{
    labels: {q_bins_json},
    datasets: [
      {{
        label: 'Raw',
        data: {raw_qhist_json},
        backgroundColor: '#93c5fd',
        borderRadius: 3,
      }},
      {{
        label: 'Scrubbed',
        data: {scrubbed_qhist_json},
        backgroundColor: '#34d399',
        borderRadius: 3,
      }},
      {{
        label: 'Downsampled',
        data: {final_qhist_json},
        backgroundColor: '#f59e0b',
        borderRadius: 3,
      }},
    ]
  }},
  options: {{
    plugins: {{ legend: {{ position: 'top' }} }},
    scales: {{ y: {{ beginAtZero: true }} }}
  }}
}});
</script>
</body>
</html>
"""

_ASSEMBLY_SECTION = """\
<h2>Assembly Statistics</h2>
<div class="grid">
  <div class="card"><div class="value">{n_contigs}</div><div class="label">Contigs</div></div>
  <div class="card"><div class="value">{total_length_kb} kb</div><div class="label">Total assembly length</div></div>
  <div class="card"><div class="value">{longest_kb} kb</div><div class="label">Longest contig</div></div>
</div>
<table>
  <thead><tr><th>Contig</th><th>Length (bp)</th><th>Status</th></tr></thead>
  <tbody>
    {contig_rows}
  </tbody>
</table>
"""


# ---------------------------------------------------------------------------
# Main report generation function
# ---------------------------------------------------------------------------

def generate_report(
    raw_reads: str,
    scrubbed_reads: str,
    final_reads: str,
    output_html: str,
    numt_stats: Optional[Dict] = None,
    downsample_stats: Optional[Dict] = None,
    assembly_path: Optional[str] = None,
    parameters: Optional[Dict] = None,
) -> None:
    """Generate the HTML report and write to output_html."""

    print("[report] Computing read statistics …")
    raw_s = read_stats(raw_reads)
    scrb_s = read_stats(scrubbed_reads)
    fin_s = read_stats(final_reads)

    numt_removed = max(0, raw_s["count"] - scrb_s["count"])
    reduction_pct = (
        round((1 - fin_s["count"] / raw_s["count"]) * 100, 1)
        if raw_s["count"] > 0 else 0
    )

    # Achieved coverage from downsampler stats
    # Use length-weighted mean across all chromosomes so a tiny chromosome
    # does not skew the displayed average.  Also build a per-contig table.
    achieved_cov = "N/A"
    per_contig_cov_rows = ""
    if downsample_stats and "achieved_coverage" in downsample_stats:
        cov_dict: Dict[str, float] = downsample_stats["achieved_coverage"]
        # Length-weighted mean requires contig lengths — use achieved bases / cov as proxy.
        # If ref_lengths are provided in stats, use them; otherwise simple mean.
        ref_lens: Optional[Dict[str, int]] = downsample_stats.get("ref_lengths")
        if ref_lens and cov_dict:
            total_len = sum(ref_lens.get(c, 0) for c in cov_dict)
            achieved_cov = (
                str(round(
                    sum(cov * ref_lens.get(c, 0) for c, cov in cov_dict.items()) / total_len, 1
                )) + "x"
            ) if total_len > 0 else "N/A"
        elif cov_dict:
            achieved_cov = str(round(sum(cov_dict.values()) / len(cov_dict), 1)) + "x"

        for chrom, cov in sorted(cov_dict.items()):
            badge = "green" if cov >= 80 else ("yellow" if cov >= 30 else "red")
            per_contig_cov_rows += (
                f'<tr><td>{chrom}</td>'
                f'<td><span class="badge badge-{badge}">{cov:.1f}x</span></td></tr>\n'
            )

    # NUMT breakdown
    numt_labels = ["Clean reads"]
    numt_values = [scrb_s["count"]]
    if numt_stats:
        for key, label in [
            ("failed_mapq", "Low MAPQ"),
            ("failed_mapped_fraction", "Poor alignment"),
            ("failed_softclip", "Chimeric (softclip)"),
            ("failed_low_coverage", "Low coverage (NUMT)"),
        ]:
            v = numt_stats.get(key, 0)
            if v > 0:
                numt_labels.append(label)
                numt_values.append(v)

    # Q-score bins
    q_bins = raw_s.get("q_bins", [f"Q{i*5}–Q{i*5+5}" for i in range(9)])

    def pad_hist(h: List[int], n: int) -> List[int]:
        return (h + [0] * n)[:n]

    n_bins = len(q_bins)
    raw_qh = pad_hist(raw_s.get("q_histogram", []), n_bins)
    scrb_qh = pad_hist(scrb_s.get("q_histogram", []), n_bins)
    fin_qh = pad_hist(fin_s.get("q_histogram", []), n_bins)

    # Assembly section
    assembly_section = ""
    if assembly_path and os.path.exists(assembly_path):
        asm = fasta_assembly_stats(assembly_path)
        contig_rows = ""
        for c in asm["contigs"]:
            badge = "green" if c["length"] > 50_000 else ("yellow" if c["length"] > 10_000 else "red")
            status = "Chromosome-scale" if c["length"] > 50_000 else ("Large" if c["length"] > 10_000 else "Fragment")
            contig_rows += (
                f'<tr><td>{c["name"]}</td>'
                f'<td>{c["length"]:,}</td>'
                f'<td><span class="badge badge-{badge}">{status}</span></td></tr>\n'
            )
        assembly_section = _ASSEMBLY_SECTION.format(
            n_contigs=asm["n_contigs"],
            total_length_kb=round(asm["total_length"] / 1000, 1),
            longest_kb=round(asm["longest"] / 1000, 1),
            contig_rows=contig_rows,
        )

    # Parameters table
    params = parameters or {}
    param_rows = "\n".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in params.items()
    )
    if not param_rows:
        param_rows = "<tr><td colspan='2'>No parameters recorded</td></tr>"

    if not per_contig_cov_rows:
        per_contig_cov_rows = "<tr><td colspan='2'>Coverage data not available</td></tr>"

    html = _HTML_TEMPLATE.format(
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        raw_count=f"{raw_s['count']:,}",
        scrubbed_count=f"{scrb_s['count']:,}",
        final_count=f"{fin_s['count']:,}",
        numt_removed=f"{numt_removed:,}",
        reduction_pct=reduction_pct,
        achieved_cov=achieved_cov,
        # table
        raw_total_mb=round(raw_s["total_bases"] / 1e6, 2),
        raw_mean_len=raw_s["mean_length"],
        raw_n50=raw_s["n50"],
        raw_mean_q=raw_s["mean_q"],
        scrubbed_total_mb=round(scrb_s["total_bases"] / 1e6, 2),
        scrubbed_mean_len=scrb_s["mean_length"],
        scrubbed_n50=scrb_s["n50"],
        scrubbed_mean_q=scrb_s["mean_q"],
        final_total_mb=round(fin_s["total_bases"] / 1e6, 2),
        final_mean_len=fin_s["mean_length"],
        final_n50=fin_s["n50"],
        final_mean_q=fin_s["mean_q"],
        # charts
        numt_labels_json=json.dumps(numt_labels),
        numt_values_json=json.dumps(numt_values),
        q_bins_json=json.dumps(q_bins),
        raw_qhist_json=json.dumps(raw_qh),
        scrubbed_qhist_json=json.dumps(scrb_qh),
        final_qhist_json=json.dumps(fin_qh),
        assembly_section=assembly_section,
        per_contig_cov_rows=per_contig_cov_rows,
        param_rows=param_rows,
    )

    with open(output_html, "w") as fh:
        fh.write(html)
    print(f"[report] HTML report written to {output_html}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate HTML QC report for MitoAssembler pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--raw-reads", required=True, dest="raw_reads")
    p.add_argument("--scrubbed-reads", required=True, dest="scrubbed_reads")
    p.add_argument("--final-reads", required=True, dest="final_reads")
    p.add_argument("--output", required=True)
    p.add_argument("--numt-stats", default=None, dest="numt_stats")
    p.add_argument("--downsample-stats", default=None, dest="downsample_stats")
    p.add_argument("--assembly", default=None)
    p.add_argument(
        "--params-json",
        default=None,
        dest="params_json",
        help="Optional path to pipeline_params.json to populate the Parameters table",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    numt_stats = None
    ds_stats = None
    parameters = None
    if args.numt_stats and os.path.exists(args.numt_stats):
        with open(args.numt_stats) as fh:
            numt_stats = json.load(fh)
    if args.downsample_stats and os.path.exists(args.downsample_stats):
        with open(args.downsample_stats) as fh:
            ds_stats = json.load(fh)
    if args.params_json and os.path.exists(args.params_json):
        with open(args.params_json) as fh:
            parameters = json.load(fh)
    generate_report(
        raw_reads=args.raw_reads,
        scrubbed_reads=args.scrubbed_reads,
        final_reads=args.final_reads,
        output_html=args.output,
        numt_stats=numt_stats,
        downsample_stats=ds_stats,
        assembly_path=args.assembly,
        parameters=parameters,
    )


if __name__ == "__main__":
    main()
