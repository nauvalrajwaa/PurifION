# MitoAssembler

Plant mitogenome assembly pipeline for ONT long reads. Performs per-chromosome filtering, downsampling, optional assembly, and mandatory graph disentanglement.

## Current Workflow

```text
raw ONT reads
      |
      v
[1] Per-chromosome PAF filter
    minimap2 (map-ont) → PAF
    Each read assigned to best-matching chromosome
    Dual-threshold: (len≥len_a AND id≥id_a) OR (len≥len_b AND id≥id_b)
      |
      v (per chromosome)
[2] Per-chromosome Downsampler
    uniform coverage normalization within each chromosome
      |
      v (per chromosome, if --run-assembler)
[3] Per-chromosome Flye assembly
    assembly.fasta + assembly_graph.gfa per chromosome
      |
      v (per chromosome, mandatory after assembly)
[4] disentangle_graph.py
    terminal-overlap graph, cycle/path search, stitched FASTA
      |
      v
[5] report.py  →  report/report.html
```

> **No scrub step.** The numt/paf scrubbing modes have been removed. The pipeline goes straight from raw reads into per-chromosome PAF filtering.

---

## Per-chromosome PAF Filter

Reads are aligned to the full reference with `minimap2 -cx map-ont`. Each read is assigned to the chromosome it aligns to best (highest score = matching\_bases × aligned\_length). A read is **kept** if it passes either threshold:

| Threshold | Length | Identity |
|-----------|--------|----------|
| A (long)  | ≥ `--paf-len-a` (5000 bp) | ≥ `--paf-id-a` (0.90) |
| B (medium)| ≥ `--paf-len-b` (3000 bp) | ≥ `--paf-id-b` (0.98) |

`identity = matching_bases / aligned_length`

Each chromosome gets its own `filtered.fastq`, which feeds its own independent downsampler and assembler.

---

## Post-Assembly: `disentangle_graph.py` (mandatory)

After every successful chromosome assembly, `disentangle_graph.py` is run automatically. It performs:

- All-vs-all oriented terminal overlap graph
- Per-contig self-terminal circularization checks
- Expected-size and any-size cycle/path search
- Stitched candidate FASTA with contig-name-based edge naming

**Edge naming:** merged path sequences are named after the actual contigs traversed:
```
contig_1__contig_2|cycle|len=490000|mean_id=89.50
```

Outputs per chromosome (under `02_assembly/<chrom>/`):
- `disentangle.terminal_overlaps.tsv`
- `disentangle.self_terminal.tsv`
- `disentangle.paths.tsv`
- `disentangle.cycles_any.tsv`
- `disentangle.path_sequences.fasta`
- `disentangle.report.txt`

---

## Requirements

### Python

```bash
pip install pysam biopython
```

### System tools

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | ≥ 2.26 | Read alignment and PAF filtering |
| samtools | ≥ 1.17 | BAM handling (downsampler) |
| flye | ≥ 2.9 | Assembly (optional) |

Recommended conda environment:

```bash
conda create -n mitoasm python=3.11 pysam biopython minimap2 samtools flye -c bioconda -c conda-forge
conda activate mitoasm
```

---

## Running with `run.sh` (Recommended)

`run.sh` handles reference fetch/prep and calls `pipeline.py`.

### A) Per-chromosome filter + downsampling only (no assembly)

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
TARGET_COV=100 \
THREADS=16 \
bash run.sh
```

### B) Full pipeline: filter + downsample + assembly + disentangle

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
RUN_ASSEMBLER=1 \
GENOME_SIZE=500k \
TARGET_COV=100 \
THREADS=16 \
bash run.sh
```

### C) Custom PAF thresholds + full pipeline

```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
PAF_LEN_A=5000 \
PAF_ID_A=0.90 \
PAF_LEN_B=3000 \
PAF_ID_B=0.98 \
RUN_ASSEMBLER=1 \
GENOME_SIZE=500k \
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

---

## Running `pipeline.py` Directly

### 1) Filter + downsampling only

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results/ \
  --threads 16 \
  --target-coverage 100
```

### 2) Full pipeline: filter + downsample + assembly + disentangle

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results/ \
  --threads 16 \
  --target-coverage 100 \
  --run-assembler \
  --assembler flye \
  --genome-size 500k
```

### 3) Custom PAF thresholds + custom disentangle params

```bash
python pipeline.py \
  --reads raw_reads.fastq \
  --reference reference/mito_ref.fasta \
  --outdir results/ \
  --threads 16 \
  --paf-len-a 5000 \
  --paf-id-a 0.90 \
  --paf-len-b 3000 \
  --paf-id-b 0.98 \
  --target-coverage 100 \
  --run-assembler \
  --assembler flye \
  --genome-size 500k \
  --disentangle-expected-size 494000 \
  --disentangle-size-tolerance 0.20 \
  --disentangle-window 500 \
  --disentangle-min-id 75 \
  --disentangle-min-ovlp 100 \
  --disentangle-max-depth 8
```

---

## Standalone Graph Analysis

### `disentangle_graph.py` on an existing GFA

```bash
python disentangle_graph.py assembly_graph.gfa \
  --window 500 \
  --min-id 75 \
  --min-ovlp 100 \
  --expected-size 494000 \
  --size-tolerance 0.20 \
  --max-depth 8 \
  --out-prefix results/disentangle
```

### Try alternate expected sizes for difficult samples

```bash
python disentangle_graph.py assembly_graph.gfa \
  --expected-size 70000 \
  --size-tolerance 0.20 \
  --out-prefix results/disentangle_70k
```

### `check_joins.py` on an existing GFA (standalone)

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

## Output Layout

```text
results/
├── 01_per_chrom/
│   ├── <chrom_1>/
│   │   ├── filtered.fastq           # PAF-filtered reads for this chromosome
│   │   ├── downsampled.fastq        # uniform-coverage downsampled reads
│   │   └── downsample_stats.json
│   └── <chrom_2>/
│       └── ...
├── 02_assembly/                     # present when --run-assembler
│   ├── <chrom_1>/
│   │   ├── assembly.fasta
│   │   ├── assembly_graph.gfa
│   │   ├── disentangle.terminal_overlaps.tsv
│   │   ├── disentangle.self_terminal.tsv
│   │   ├── disentangle.paths.tsv
│   │   ├── disentangle.cycles_any.tsv
│   │   ├── disentangle.path_sequences.fasta
│   │   └── disentangle.report.txt
│   └── <chrom_2>/
│       └── ...
├── report/
│   └── report.html
└── pipeline_params.json
```

---

## Key CLI Parameters

### PAF Filter

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--paf-len-a` | `5000` | Threshold A: minimum alignment length |
| `--paf-id-a` | `0.90` | Threshold A: minimum identity |
| `--paf-len-b` | `3000` | Threshold B: minimum alignment length |
| `--paf-id-b` | `0.98` | Threshold B: minimum identity |

### Downsampler

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--target-coverage` | `100` | Target normalized coverage per chromosome |
| `--bin-size` | `500` | Bin size in bp for greedy selection |
| `--spanning-fraction` | `0.80` | Always retain reads spanning ≥ this fraction |

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run-assembler` | off | Enable per-chromosome Flye assembly |
| `--assembler` | `flye` | Assembler to use |
| `--genome-size` | `500k` | Flye genome size estimate |

### Disentangle (auto-applied after assembly)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--disentangle-window` | `500` | Terminal window size (bp) |
| `--disentangle-min-id` | `75.0` | Minimum terminal overlap identity (%) |
| `--disentangle-min-ovlp` | `100` | Minimum overlap length (bp) |
| `--disentangle-top-out` | `5` | Top N edges per source+orientation |
| `--disentangle-expected-size` | `494000` | Expected circular chromosome size |
| `--disentangle-size-tolerance` | `0.20` | Size tolerance fraction (±20%) |
| `--disentangle-max-depth` | `8` | Max DFS depth for path search |
| `--disentangle-max-cycles` | `30` | Max expected-size cycles to report |
| `--disentangle-max-linear` | `50` | Max linear paths to report |
| `--disentangle-max-cycles-any` | `200` | Max any-size cycles to report |
| `--disentangle-fasta-top` | `20` | Top N sequences to write to FASTA |

---

## Notes for Difficult Samples

- If no full circular path is found, inspect `disentangle.cycles_any.tsv` and `disentangle.paths.tsv` for partial loops or strong local joins.
- Try multiple expected sizes with `--disentangle-expected-size` to account for size variation.
- Per-chromosome `filtered.fastq` and `downsampled.fastq` can be reused for manual assembly attempts.
- Edge names in `disentangle.path_sequences.fasta` encode the merged contigs: `contig_A__contig_B|kind|len=...|mean_id=...`.
