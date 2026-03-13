# MitoAssembler

Plant mitogenome assembly pipeline for ONT long reads (Scheme 4).  
Multi-chromosome PAF filter → per-chromosome downsampling → per-chromosome assembly → mandatory disentangle graph → HTML report.

---

## Why Scheme 4?

Plant mitogenomes have two properties that break simpler schemes:

1. **Inter-chromosomal repeat sharing** — the same repeat unit (tandem, inverted, dispersed) appears on multiple chromosomes. A read spanning a repeat-to-unique boundary is genuinely relevant to *both* chromosome assemblies.
2. **Uneven coverage** — AT-rich regions, structural variant sites, and chromosome junctions have naturally low local depth.

| Scheme | Problem |
|--------|---------|
| Scheme 1: per-chrom filter (single best assignment) | Reads near shared repeats get hard-assigned to one chrom; the other assembly never sees them → gaps at repeat boundaries |
| Scheme 2: universal filter → per-chrom filter | Double filtering = double opportunity to drop real mito reads |
| Scheme 3: per-chrom filter → merge → per-chrom assemble | Logically contradictory; chromosome assignments are lost after merge |
| **Scheme 4 (this pipeline)** | Reads assigned to **every** chromosome they validly align to; inter-chromosomal repeat-spanning reads reach all relevant assemblies |

---

## Workflow

```text
raw ONT reads
      |
      v
[1] Multi-chromosome PAF filter   (minimap2 --secondary=yes -N 10)
      |
      |  Each read is placed into EVERY chromosome pool where it
      |  passes the dual-threshold filter (identity + length).
      |  A single read can appear in multiple pools.
      |
      v
[2] Per-chromosome Downsampler    (spanning reads prioritised first,
      |                            then greedy bin fill to --target-coverage)
      |                            → skip entirely with --skip-downsample
      |                              to pass full-depth reads to assembler
      |                              (mirrors empirically successful direct-align method)
      v
[3] Per-chromosome Flye assembly  (optional, --run-assembler)
      |
      v
[4] Per-chromosome Disentangle    (MANDATORY when assembler runs)
      |                            disentangle_graph.py — terminal-overlap
      |                            graph, cycle/path search, stitched FASTA
      |                            Edges named after merged contigs:
      |                            contig_A__contig_B|cycle|len=...|mean_id=...
      v
[5] HTML Report                   (report.py)
```

---

## Multi-Chromosome PAF Filter — How It Works

```
minimap2 -cx map-ont --secondary=yes -N 10 <reference.fasta> <reads.fastq>
                 ↓
         PAF output (primary + secondary alignments)
                 ↓
   For each alignment row:
     identity = matches / aligned_length
     passes = (aligned_length >= len_a  AND  identity >= id_a)
           OR (aligned_length >= len_b  AND  identity >= id_b)
                 ↓
     IF passes → add read_name to chrom_pool[target_chrom]
     (same read can be added to MULTIPLE chrom pools)
                 ↓
   Write per-chrom filtered.fastq from the union of all needed reads
```

**Key parameter `--paf-no-secondary`:** disables secondary alignments, reverting to single-best-chrom assignment (Scheme 1 behaviour). Use only for benchmarking/comparison.

---

## PAF Threshold Grid — Recommended Settings

Based on empirical results (direct alignment → assembly without downsampling).
Profiles are ordered from **strictest** (top) to **most relaxed** (bottom).

| # | Profile | `--paf-len-a` | `--paf-id-a` | `--paf-len-b` | `--paf-id-b` | When to use |
|---|---------|--------------|-------------|--------------|-------------|-------------|
| 1 | **Ultra-strict** | `8000` | `0.92` | `5000` | `0.98` | Very high depth (>500×), heavy NUMT contamination, only the cleanest mito reads |
| 2 | **Strict (default)** | `5000` | `0.90` | `3000` | `0.98` | Standard samples, clean sequencing — empirically 2+ chrom success |
| 3 | **Strict-B relaxed** | `5000` | `0.90` | `2000` | `0.95` | Standard samples, slightly softer short-read threshold |
| 4 | **Moderate** | `4000` | `0.85` | `2000` | `0.92` | Moderate quality ONT runs, slightly degraded DNA |
| 5 | **Relaxed-A** | `5000` | `0.80` | `1500` | `0.90` | Mixed-quality libraries, moderate NUMT risk |
| 6 | **Relaxed (empirical fallback)** | `5000` | `0.70` | `1000` | `0.85` | Difficult samples — empirically 1-chrom success |
| 7 | **Permissive** | `3000` | `0.80` | `1000` | `0.90` | Low-depth samples, short-read contamination |
| 8 | **Permissive-relaxed** | `3000` | `0.75` | `800` | `0.85` | Very low depth, degraded DNA, last resort before giving up on identity filter |
| 9 | **Length-biased** | `10000` | `0.70` | `5000` | `0.80` | Very long reads (N50 > 20 kb), prioritise spanning reads over identity |
| 10 | **Ultra-permissive** | `2000` | `0.70` | `500` | `0.80` | Severely degraded samples, old herbarium specimens, minimal quality requirement |
| 11 | **Short-read friendly** | `1000` | `0.85` | `500` | `0.90` | Short ONT reads (N50 < 5 kb) or hybrid short+long datasets |
| 12 | **Identity-first** | `3000` | `0.95` | `1000` | `0.98` | Highly accurate reads (Dorado super-accuracy model), low NUMT risk |

> **How to choose:**
> 1. Start with **Strict (profile 2)** — this is the empirically tested default.
> 2. If assembly has gaps or missing chromosomes → step down to **Moderate (4)** or **Relaxed-A (5)**.
> 3. If the assembly graph is tangled with NUMT contigs → step up to **Ultra-strict (1)**.
> 4. For difficult/degraded samples where strict fails → use **Relaxed fallback (6)** which had empirical 1-chrom success.
> 5. For very short reads or old specimens → try **Short-read friendly (11)** or **Ultra-permissive (10)**.
> 6. If using Dorado super-accuracy basecalling → **Identity-first (12)** captures only high-confidence alignments.

### Downsampler strategy

| Coverage situation | Recommendation |
|-------------------|----------------|
| Raw ≤ 200× | Use `--skip-downsample` — no reads to spare |
| Raw 200–500× | Default (`--target-coverage 200`) is safe |
| Raw > 500× | Default is fine; raise to `--target-coverage 300` if assembly has gaps |
| Empirically tested method | `--skip-downsample` — passes full-depth filtered reads, matches direct-align success |

---

## Output Layout

```text
results/
├── 01_per_chrom/
│   ├── <chrom_1>/
│   │   ├── filtered.fastq          # reads passing PAF filter for this chrom
│   │   ├── downsampled.fastq       # coverage-normalised subset (or copy of filtered if --skip-downsample)
│   │   └── downsample_stats.json
│   ├── <chrom_2>/
│   │   └── ...
│   └── ...
├── 02_assembly/
│   ├── <chrom_1>/
│   │   ├── assembly.fasta
│   │   ├── assembly_graph.gfa
│   │   ├── disentangle.report.txt
│   │   ├── disentangle.terminal_overlaps.tsv
│   │   ├── disentangle.paths.tsv
│   │   ├── disentangle.cycles_any.tsv
│   │   ├── disentangle.self_terminal.tsv
│   │   └── disentangle.path_sequences.fasta
│   ├── <chrom_2>/
│   │   └── ...
│   └── ...
├── report/
│   └── report.html
└── pipeline_params.json
```

---

## Requirements

### Python

```bash
pip install pysam biopython
```

### System tools

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | ≥ 2.26 | Read alignment + PAF |
| samtools | ≥ 1.17 | BAM handling (downsampler) |
| flye | ≥ 2.9 | Assembly (optional) |

Recommended conda environment:

```bash
conda create -n mitoasm python=3.11 pysam biopython minimap2 samtools flye \
  -c bioconda -c conda-forge
conda activate mitoasm
```

---

## Quick Start Examples

### A) Filter + downsample only (no assembly)

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results/ \
  --threads 16
```

### B) Full pipeline — strict thresholds + downsample + assembly

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results/ \
  --threads 16 \
  --paf-len-a 5000 --paf-id-a 0.90 \
  --paf-len-b 3000 --paf-id-b 0.98 \
  --target-coverage 200 \
  --run-assembler \
  --genome-size 500k
```

### C) Mirror empirically successful method — skip downsampler, strict thresholds

This replicates the direct-alignment → full-depth assembly approach that gave 2 successful samples:

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_strict_full/ \
  --threads 16 \
  --paf-len-a 5000 --paf-id-a 0.90 \
  --paf-len-b 3000 --paf-id-b 0.98 \
  --skip-downsample \
  --run-assembler \
  --genome-size 500k
```

### D) Relaxed thresholds + skip downsampler (1-sample fallback)

Mirrors the `L5000_id0.70_L1000_id0.85` profile that rescued the difficult sample:

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_relaxed_full/ \
  --threads 16 \
  --paf-len-a 5000 --paf-id-a 0.70 \
  --paf-len-b 1000 --paf-id-b 0.85 \
  --skip-downsample \
  --run-assembler \
  --genome-size 500k
```

### E) With downsampling at high coverage (>300× raw)

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_downsampled/ \
  --threads 16 \
  --paf-len-a 5000 --paf-id-a 0.90 \
  --paf-len-b 3000 --paf-id-b 0.98 \
  --target-coverage 200 \
  --spanning-fraction 0.80 \
  --run-assembler \
  --genome-size 500k
```

### F) Single-best-chrom mode (Scheme 1 comparison / benchmarking)

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_scheme1/ \
  --threads 16 \
  --paf-no-secondary \
  --skip-downsample \
  --run-assembler
```

### G) Custom disentangle parameters (large genome / difficult graph)

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_large/ \
  --threads 16 \
  --skip-downsample \
  --run-assembler \
  --genome-size 800k \
  --disentangle-expected-size 700000 \
  --disentangle-size-tolerance 0.25 \
  --disentangle-max-depth 12 \
  --disentangle-max-cycles 50 \
  --disentangle-fasta-top 30
```

---

## Running with `run.sh` (Recommended)

`run.sh` handles reference fetch/prep and calls `pipeline.py`.

### A) Local reference, full pipeline, skip downsampler (empirically tested)

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
RUN_ASSEMBLER=1 \
GENOME_SIZE=500k \
SKIP_DOWNSAMPLE=1 \
THREADS=16 \
bash run.sh
```

### B) Local reference, full pipeline, with downsampling

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
RUN_ASSEMBLER=1 \
GENOME_SIZE=500k \
TARGET_COV=200 \
THREADS=16 \
bash run.sh
```

### C) NCBI reference fetch

```bash
READS=raw_reads.fastq \
REF_SOURCE=ncbi \
THREADS=16 \
bash run.sh
```

Set `ACCESSIONS=(...)` inside `run.sh` for your species/chromosomes.

---

## Full CLI Reference

```
python pipeline.py --help
```

### Required

| Parameter | Description |
|-----------|-------------|
| `--reads` | Input reads (FASTQ or FASTA, ONT) |
| `--reference` | Multi-chromosome mitogenome reference FASTA |
| `--outdir` | Output directory (created if absent) |

### General

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | `4` | CPU threads for minimap2 / Flye |
| `--report-name` | `report.html` | HTML report filename inside `report/` |

### PAF Filter (Scheme 4)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--paf-len-a` | `5000` | Threshold A: minimum aligned length (bp) |
| `--paf-id-a` | `0.90` | Threshold A: minimum identity (0–1) |
| `--paf-len-b` | `3000` | Threshold B: minimum aligned length (bp) |
| `--paf-id-b` | `0.98` | Threshold B: minimum identity (0–1) |
| `--paf-no-secondary` | *(off)* | Disable secondary alignments → single-best-chrom assignment (Scheme 1) |

A read passes if: `(aln_len ≥ len_a AND identity ≥ id_a)` **OR** `(aln_len ≥ len_b AND identity ≥ id_b)`

### Downsampler

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--target-coverage` | `200` | Target normalised coverage per chromosome |
| `--bin-size` | `500` | Bin size for greedy coverage fill (bp) |
| `--spanning-fraction` | `0.80` | Reads spanning ≥ this fraction of chrom length are always kept first |
| `--skip-downsample` | *(off)* | Skip downsampler entirely — full-depth filtered reads go straight to assembler. Mirrors empirically successful direct-align method. |

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run-assembler` | *(off)* | Enable per-chromosome Flye assembly |
| `--assembler` | `flye` | Assembler to use (currently: `flye`) |
| `--genome-size` | `500k` | Flye genome size estimate per chromosome |

### Disentangle Graph (auto-runs with `--run-assembler`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--disentangle-window` | `500` | Terminal window size (bp) for overlap search |
| `--disentangle-min-id` | `75.0` | Minimum terminal overlap identity (%) |
| `--disentangle-min-ovlp` | `100` | Minimum terminal overlap length (bp) |
| `--disentangle-top-out` | `5` | Max edges retained per source+orientation |
| `--disentangle-expected-size` | `494000` | Expected circular chromosome size (bp) |
| `--disentangle-size-tolerance` | `0.20` | ± size tolerance fraction (e.g. 0.20 = ±20%) |
| `--disentangle-max-depth` | `8` | Max DFS depth for path/cycle search |
| `--disentangle-max-cycles` | `30` | Max cycles within expected size to report |
| `--disentangle-max-linear` | `50` | Max linear paths within expected size to report |
| `--disentangle-max-cycles-any` | `200` | Max any-size cycles to report |
| `--disentangle-fasta-top` | `20` | Top N stitched path sequences to write to FASTA |

---

## Standalone Tools

### `disentangle_graph.py` — Deep graph analysis on any GFA

```bash
python disentangle_graph.py assembly_graph.gfa \
  --window 500 \
  --min-id 75 \
  --min-ovlp 100 \
  --expected-size 494000 \
  --size-tolerance 0.20 \
  --max-depth 8 \
  --out-prefix results/disentangle_full
```

Outputs:
- `PREFIX.terminal_overlaps.tsv` — all-vs-all oriented terminal overlaps
- `PREFIX.self_terminal.tsv` — per-contig self-circularization candidates
- `PREFIX.paths.tsv` — linear paths within expected size
- `PREFIX.cycles_any.tsv` — cycles of any size
- `PREFIX.path_sequences.fasta` — stitched candidate sequences
  - Named: `contig_A__contig_B|cycle|len=...|mean_id=...`
- `PREFIX.report.txt` — consolidated text summary

Try alternate expected sizes for difficult samples:

```bash
python disentangle_graph.py assembly_graph.gfa \
  --expected-size 70000 \
  --size-tolerance 0.30 \
  --out-prefix results/disentangle_70k
```

### `check_joins.py` — Quick terminal join scan on any GFA

```bash
python check_joins.py assembly_graph.gfa \
  --window 500 \
  --min-id 75 \
  --min-ovlp 100 \
  --expected-size 494000 \
  --size-tolerance 0.20 \
  --out candidate_joins.tsv
```

---

## Troubleshooting Difficult Samples

| Problem | Action |
|---------|--------|
| No circular path found | Check `disentangle.cycles_any.tsv`; lower `--disentangle-min-id`; try multiple `--disentangle-expected-size` values |
| Assembly gaps at repeat boundaries | Ensure `--paf-no-secondary` is NOT set; try `--skip-downsample` to keep full-depth reads |
| Low read count per chromosome | Switch to **Relaxed** PAF grid (`--paf-len-a 5000 --paf-id-a 0.70 --paf-len-b 1000 --paf-id-b 0.85`) |
| Very slow assembly | Raise `--paf-len-a` / `--paf-id-a` or lower `--target-coverage`; do NOT use `--skip-downsample` at >500× raw |
| NUMT contamination in graph | Use **Ultra-strict** PAF grid (`--paf-len-a 8000 --paf-id-a 0.92 --paf-len-b 5000 --paf-id-b 0.98`) |
| Downsampler drops too many reads | Use `--skip-downsample`; or raise `--spanning-fraction` to 0.90 |
