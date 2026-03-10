# MitoAssembler

Plant mitogenome assembly pipeline for ONT long reads, with two filtering modes and optional post-assembly graph-join analysis.

## Current Workflow

```text
raw ONT reads
      |
      |  choose one filtering mode
      |  - numt: scrubber-based filtering
      |  - paf : dual-threshold PAF filtering
      v
[1] Filtered reads
      |
      v
[2] downsampler.py (uniform coverage normalization)
      |
      v
[3] (optional) Flye assembly
      |
      +--> (optional) check_joins.py on assembly_graph.gfa
      |
      +--> (optional, standalone) disentangle_graph.py for deep terminal-overlap/cycle search
      v
[4] report.py (HTML report)
```

## Two Filter Modes (Pipeline)

The pipeline supports two filtering strategies via `--filter-mode` (or `FILTER_MODE` in `run.sh`):

1. `numt`
- Uses `numt_scrubber.py` logic (MAPQ/mapped-fraction/soft-clip/coverage-based filtering).
- Recommended default for most datasets.

2. `paf`
- Uses old dual-threshold PAF filter strategy (compatible with `old_methods` behavior).
- Tunable with:
  - `--paf-len-a` / `--paf-id-a`
  - `--paf-len-b` / `--paf-id-b`

## Post-Assembly Analysis

### 1) `check_joins.py` (integrated into pipeline)

When assembly is enabled, you can automatically run terminal join analysis on Flye graph:

- Enable with `--run-check-joins` (or `RUN_CHECK_JOINS=1` in `run.sh`).
- Input: `03_assembly/assembly_graph.gfa`
- Outputs:
  - `03_assembly/candidate_joins.tsv`
  - `03_assembly/check_joins_summary.txt`

### 2) `disentangle_graph.py` (standalone, deeper analysis)

Use after assembly (or on any GFA) for:
- all-vs-all oriented terminal overlap graph,
- per-contig self-terminal circularization checks,
- expected-size and any-size cycle/path search,
- stitched candidate FASTA outputs,
- consolidated text report.

Typical outputs (`--out-prefix PREFIX`):
- `PREFIX.terminal_overlaps.tsv`
- `PREFIX.self_terminal.tsv`
- `PREFIX.paths.tsv`
- `PREFIX.cycles_any.tsv`
- `PREFIX.path_sequences.fasta`
- `PREFIX.report.txt`

## Requirements

### Python

```bash
pip install pysam biopython
```

### System tools

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | >= 2.26 | Read alignment |
| samtools | >= 1.17 | BAM handling |
| flye | >= 2.9 | Assembly (optional) |

Recommended conda environment:

```bash
conda create -n mitoasm python=3.11 pysam biopython minimap2 samtools flye -c bioconda -c conda-forge
conda activate mitoasm
```

## Running with `run.sh` (Recommended)

`run.sh` handles reference fetch/prep and calls `pipeline.py`.

### A) Default NUMT mode + downsampling (no assembly)

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
FILTER_MODE=numt \
TARGET_COV=100 \
THREADS=16 \
bash run.sh
```

### B) PAF mode (old_methods style) + downsampling

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
FILTER_MODE=paf \
PAF_LEN_A=5000 \
PAF_ID_A=0.90 \
PAF_LEN_B=3000 \
PAF_ID_B=0.98 \
TARGET_COV=100 \
THREADS=16 \
bash run.sh
```

### C) Full pipeline: PAF mode + assembly + check_joins

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
FILTER_MODE=paf \
RUN_ASSEMBLER=1 \
GENOME_SIZE=500k \
RUN_CHECK_JOINS=1 \
JOINS_WINDOW=500 \
JOINS_MIN_ID=75 \
JOINS_MIN_OVLP=100 \
JOINS_EXPECTED_SIZE=494000 \
JOINS_SIZE_TOLERANCE=0.20 \
TARGET_COV=100 \
THREADS=16 \
bash run.sh
```

### D) NCBI reference fetch mode

```bash
READS=raw_reads.fastq \
REF_SOURCE=ncbi \
THREADS=16 \
bash run.sh
```

Set `ACCESSIONS=(...)` inside `run.sh` for your species/chromosomes.

## Running `pipeline.py` Directly

### 1) NUMT mode + downsampling

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_numt \
  --threads 16 \
  --filter-mode numt \
  --target-coverage 100
```

### 2) PAF mode + downsampling

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_paf \
  --threads 16 \
  --filter-mode paf \
  --paf-len-a 5000 \
  --paf-id-a 0.90 \
  --paf-len-b 3000 \
  --paf-id-b 0.98 \
  --target-coverage 100
```

### 3) Full run + assembly + integrated check_joins

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results_full \
  --threads 16 \
  --filter-mode paf \
  --paf-len-a 5000 \
  --paf-id-a 0.90 \
  --paf-len-b 3000 \
  --paf-id-b 0.98 \
  --target-coverage 100 \
  --run-assembler \
  --assembler flye \
  --genome-size 500k \
  --run-check-joins \
  --joins-window 500 \
  --joins-min-id 75 \
  --joins-min-ovlp 100 \
  --joins-expected-size 494000 \
  --joins-size-tolerance 0.20
```

## Standalone Graph Analysis Commands

### A) check_joins on an existing GFA

```bash
python check_joins.py old_methods/assembly_graph.gfa \
  --window 500 \
  --min-id 75 \
  --min-ovlp 100 \
  --expected-size 494000 \
  --size-tolerance 0.20 \
  --out old_methods/candidate_joins.tsv
```

### B) disentangle_graph on an existing GFA

```bash
python disentangle_graph.py old_methods/assembly_graph.gfa \
  --window 500 \
  --min-id 75 \
  --min-ovlp 100 \
  --expected-size 494000 \
  --size-tolerance 0.20 \
  --max-depth 8 \
  --out-prefix old_methods/disentangle_full
```

### C) Try alternate expected sizes for difficult samples

```bash
python disentangle_graph.py old_methods/assembly_graph.gfa \
  --expected-size 70000 \
  --size-tolerance 0.20 \
  --out-prefix old_methods/disentangle_70k
```

## Output Layout

```text
results/
+-- 01_numt_scrubbed/            # filtered reads output directory
|   +-- scrubbed.fastq
|   +-- numt_stats.json
+-- 02_downsampled/
|   +-- downsampled.fastq
|   +-- downsample_stats.json
+-- 03_assembly/                 # present when --run-assembler
|   +-- assembly.fasta
|   +-- assembly_graph.gfa
|   +-- candidate_joins.tsv      # when --run-check-joins
|   +-- check_joins_summary.txt  # when --run-check-joins
+-- report/
|   +-- report.html
+-- pipeline_params.json
```

`disentangle_graph.py` outputs are written wherever `--out-prefix` points.

## Key CLI Parameters

| Group | Parameter | Default | Description |
|------|-----------|---------|-------------|
| filter | `--filter-mode` | `numt` | Filtering strategy: `numt` or `paf` |
| paf | `--paf-len-a` | `5000` | PAF threshold A length |
| paf | `--paf-id-a` | `0.90` | PAF threshold A identity |
| paf | `--paf-len-b` | `3000` | PAF threshold B length |
| paf | `--paf-id-b` | `0.98` | PAF threshold B identity |
| numt | `--min-mapq` | `20` | Minimum MAPQ |
| numt | `--min-mapped-fraction` | `0.80` | Minimum mapped fraction |
| numt | `--max-softclip-fraction` | `0.20` | Maximum softclip fraction |
| numt | `--coverage-window` | `500` | Coverage window size |
| numt | `--coverage-threshold-ratio` | `0.20` | Low-coverage threshold ratio |
| downsample | `--target-coverage` | `100` | Target normalized coverage |
| downsample | `--bin-size` | `500` | Bin size (bp) |
| downsample | `--spanning-fraction` | `0.80` | Always keep spanning reads above this fraction |
| assembly | `--run-assembler` | off | Enable Flye |
| assembly | `--genome-size` | `500k` | Flye genome size estimate |
| joins | `--run-check-joins` | off | Run `check_joins.py` after assembly |
| joins | `--joins-window` | `500` | Join terminal window |
| joins | `--joins-min-id` | `75.0` | Min join identity (%) |
| joins | `--joins-min-ovlp` | `100` | Min join overlap (bp) |
| joins | `--joins-expected-size` | `494000` | Expected circular size |
| joins | `--joins-size-tolerance` | `0.20` | Size tolerance fraction |

## Notes for Difficult Samples

- If full circular mt path is not found, inspect `candidate_joins.tsv` and `disentangle` outputs for strong local loops/joins.
- Try multiple expected sizes with `disentangle_graph.py`.
- Keep both filter modes in your comparison matrix (`numt` vs `paf`) before deciding final assembly strategy.
