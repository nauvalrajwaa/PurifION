# MitoAssembler

Plant mitogenome assembly pipeline for ONT long reads.

## Pipeline overview

```
raw ONT reads
      │
      ▼
[1] numt_scrubber.py   ← removes NUMTs and chimeric reads
      │
      ▼
[2] downsampler.py     ← normalises coverage, keeps longest/highest-Q reads
      │
      ▼
[3] (optional) Flye assembly
      │
      ▼
[4] report.py          ← self-contained HTML QC report
```

## Requirements

### Python
```
pip install pysam
```

### System tools (install via conda or system package manager)
| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | ≥ 2.26 | Read alignment |
| samtools | ≥ 1.17 | BAM handling |
| flye | ≥ 2.9 | Assembly (optional) |

Recommended install via conda:
```bash
conda create -n mitoasm python=3.11 pysam minimap2 samtools flye -c bioconda -c conda-forge
conda activate mitoasm
```

## Quick start

### Using run.sh (recommended)

`run.sh` handles reference fetching and then calls `pipeline.py`.
All variables can be overridden inline without editing the file.

**NCBI mode — download chromosomes by accession number:**
```bash
# Edit the ACCESSIONS array in run.sh, then:
READS=raw_reads.fastq \
THREADS=16 \
TARGET_COV=100 \
bash run.sh
```
Or override accessions inline (space-separated, quoted as a single string):
```bash
READS=raw_reads.fastq \
REF_SOURCE=ncbi \
ACCESSIONS=(NC_002511.2 NC_002982.1 NC_003997.1) \
THREADS=16 \
bash run.sh
```

**Local mode — use a FASTA already on disk:**
```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
THREADS=16 \
TARGET_COV=100 \
bash run.sh
```

**Cat mode — concatenate multiple per-chromosome FASTAs:**
```bash
READS=raw_reads.fastq \
REF_SOURCE=cat \
REF_FILES="chr1.fasta chr2.fasta chr3.fasta" \
THREADS=16 \
bash run.sh
```

**With assembly (Flye):**
```bash
READS=raw_reads.fastq \
REF_SOURCE=local \
REF_PATH=reference/mito_ref.fasta \
RUN_ASSEMBLER=1 \
GENOME_SIZE=500k \
THREADS=16 \
bash run.sh
```

### Full pipeline without run.sh
```bash
python pipeline.py \
    --reads raw_reads.fastq \
    --reference mito_ref.fasta \
    --outdir results/ \
    --threads 16 \
    --target-coverage 100 \
    --run-assembler \
    --genome-size 500k
```

### Individual tools

**NUMT Scrubber only:**
```bash
python numt_scrubber.py \
    --reads raw_reads.fastq \
    --reference mito_ref.fasta \
    --output scrubbed.fastq \
    --threads 8 \
    --stats numt_stats.json
```

**Downsampler only:**
```bash
python downsampler.py \
    --reads scrubbed.fastq \
    --reference mito_ref.fasta \
    --output downsampled.fastq \
    --target-coverage 100 \
    --threads 8 \
    --stats downsample_stats.json
```

**Report only:**
```bash
python report.py \
    --raw-reads raw_reads.fastq \
    --scrubbed-reads scrubbed.fastq \
    --final-reads downsampled.fastq \
    --numt-stats numt_stats.json \
    --downsample-stats downsample_stats.json \
    --params-json pipeline_params.json \
    --assembly assembly.fasta \
    --output report.html
```

## Output layout

```
results/
├── 01_numt_scrubbed/
│   ├── scrubbed.fastq       # reads after NUMT removal
│   └── numt_stats.json      # per-filter breakdown
├── 02_downsampled/
│   ├── downsampled.fastq    # assembly-ready reads (~100x uniform coverage)
│   └── downsample_stats.json
├── 03_assembly/             # only when --run-assembler
│   └── assembly.fasta
├── report/
│   └── report.html          # self-contained QC report
└── pipeline_params.json     # full parameter record
```

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | 4 | CPU threads |
| `--min-mapq` | 20 | Minimum mapping quality for NUMT filter |
| `--min-mapped-fraction` | 0.80 | Minimum fraction of read that must align |
| `--max-softclip-fraction` | 0.20 | Maximum soft-clip fraction (chimeric indicator) |
| `--coverage-window` | 500 | Sliding window size (bp) for coverage estimation |
| `--coverage-threshold-ratio` | 0.20 | Fraction of baseline below which a region is NUMT |
| `--target-coverage` | 100 | Target mean coverage for downsampler |
| `--bin-size` | 500 | Reference bin size (bp) for coverage normalisation |
| `--spanning-fraction` | 0.80 | Reads covering ≥ this fraction of genome are always kept |

## Algorithm details

### NUMT Scrubber
Reads failing **any** of four filters are removed:
1. **MAPQ < min-mapq** — poorly anchored alignments
2. **Mapped fraction < min-mapped-fraction** — reads with large unaligned tails (chimeric)
3. **Soft-clip fraction > max-softclip-fraction** — reads with large nuclear flanking sequences
4. **Local coverage < baseline × threshold-ratio** — reads mapping to low-coverage regions typical of NUMT insertion sites

### Uniform Downsampler
1. Maps reads to the mitochondrial reference.
2. Divides the genome into bins of `--bin-size` bp.
3. For each bin, ranks candidate reads by **(length DESC, Q-score DESC)** and greedily fills until `target_coverage × bin_size` bases are accumulated.
4. Reads spanning ≥ `spanning-fraction` of the genome are **always retained** to ensure circular topology resolution.
