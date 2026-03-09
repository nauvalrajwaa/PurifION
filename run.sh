#!/usr/bin/env bash
# =============================================================================
# run.sh — MitoAssembler: fetch reference + run pipeline
# =============================================================================
# Usage:
#   bash run.sh                        # uses defaults below
#   READS=my.fastq bash run.sh         # override any variable inline
#
# Multi-chromosome plant mitogenomes are handled automatically.
# The --reference flag accepts a single FASTA that may contain ANY number of
# sequences (chromosomes).  Each chromosome gets its own:
#   - coverage baseline (NUMT filter)
#   - independent bin-filling budget (downsampler)
#   - spanning-read detection
#   - per-chromosome coverage row in the HTML report
#
# Reference fetch strategy (pick one by setting REF_SOURCE below):
#   ncbi   — download from NCBI using accession numbers (requires datasets CLI
#             or efetch; falls back to efetch if datasets is not found)
#   local  — skip download, use a FASTA already on disk (set REF_PATH)
#   cat    — concatenate multiple local FASTAs into one multi-chromosome FASTA
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# ① USER CONFIGURATION — edit this section
# ---------------------------------------------------------------------------

# Input ONT reads (FASTQ or FASTA, may be gzipped)
READS="${READS:-raw_reads.fastq}"

# Reference source: "ncbi" | "local" | "cat"
REF_SOURCE="${REF_SOURCE:-ncbi}"

# --- ncbi mode ---
# NCBI accession numbers for each mitochondrial chromosome.
# Plant mitogenomes often have 2–7 chromosomes published as separate records.
# Example below uses Beta vulgaris (sugar beet) mitogenome: 3 chromosomes.
# Replace with your species' accessions from:
#   https://www.ncbi.nlm.nih.gov/genome/organelle/
ACCESSIONS=(
    "NC_002511.2"   # chromosome 1  (main)
    "NC_002982.1"   # chromosome 2
    "NC_003997.1"   # chromosome 3
    # add more as needed
)

# --- local mode ---
# Path to an existing multi-chromosome FASTA (already merged, if needed)
REF_PATH="${REF_PATH:-reference/mito_ref.fasta}"

# --- cat mode ---
# Space-separated list of per-chromosome FASTA files to concatenate
REF_FILES="${REF_FILES:-chr1.fasta chr2.fasta chr3.fasta}"

# Output directory
OUTDIR="${OUTDIR:-results}"

# Number of CPU threads
THREADS="${THREADS:-8}"

# Assembler: set to 1 to run Flye after downsampling, 0 to skip
RUN_ASSEMBLER="${RUN_ASSEMBLER:-0}"

# Estimated total mitogenome size (used by Flye; ignored when RUN_ASSEMBLER=0)
# Typical plant mitogenomes: 200k–2m
GENOME_SIZE="${GENOME_SIZE:-500k}"

# Target coverage for downsampler
TARGET_COV="${TARGET_COV:-100}"

# ---------------------------------------------------------------------------
# ② PATHS & HELPERS
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${OUTDIR}/reference"
MERGED_REF="${REF_DIR}/mito_combined.fasta"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

check_tool() {
    command -v "$1" &>/dev/null || die "'$1' not found. Install with: $2"
}

# ---------------------------------------------------------------------------
# ③ PREFLIGHT CHECKS
# ---------------------------------------------------------------------------

log "Checking required tools …"
check_tool python   "conda install python"
check_tool minimap2 "conda install -c bioconda minimap2"
check_tool samtools "conda install -c bioconda samtools"
[ -f "${SCRIPT_DIR}/pipeline.py" ] || die "pipeline.py not found in ${SCRIPT_DIR}"
[ -f "${READS}" ]                  || die "Reads file not found: ${READS}"

if [ "${RUN_ASSEMBLER}" = "1" ]; then
    check_tool flye "conda install -c bioconda flye"
fi

python -c "import pysam" 2>/dev/null \
    || die "pysam not installed. Run: pip install pysam  (or conda install -c bioconda pysam)"

# ---------------------------------------------------------------------------
# ④ FETCH / PREPARE REFERENCE
# ---------------------------------------------------------------------------

mkdir -p "${REF_DIR}"

case "${REF_SOURCE}" in

# ── NCBI: download each accession and concatenate ──────────────────────────
ncbi)
    log "Fetching reference sequences from NCBI …"
    log "Accessions: ${ACCESSIONS[*]}"

    FETCHED_FILES=()

    for ACC in "${ACCESSIONS[@]}"; do
        OUT_FASTA="${REF_DIR}/${ACC}.fasta"

        if [ -f "${OUT_FASTA}" ]; then
            log "  ${ACC} — already downloaded, skipping"
        else
            log "  Downloading ${ACC} …"

            # Try NCBI datasets CLI first (modern, faster).
            # NOTE: 'datasets download genome' works for full assembly accessions
            # (GCA_/GCF_).  For individual nucleotide records (NC_/MK_/etc.) it
            # will silently produce an empty file; the efetch/curl fallback below
            # handles those correctly.
            if command -v datasets &>/dev/null; then
                datasets download genome accession "${ACC}" \
                    --filename "${REF_DIR}/${ACC}.zip" \
                    --no-progressbar 2>/dev/null \
                && unzip -o "${REF_DIR}/${ACC}.zip" \
                         -d "${REF_DIR}/${ACC}_dl" >/dev/null \
                && find "${REF_DIR}/${ACC}_dl" -name "*.fna" \
                         -exec cat {} \; > "${OUT_FASTA}" \
                && rm -rf "${REF_DIR}/${ACC}.zip" "${REF_DIR}/${ACC}_dl" \
                || true
            fi

            # Fallback: efetch (NCBI E-utilities, always available via conda)
            if [ ! -s "${OUT_FASTA}" ]; then
                if command -v efetch &>/dev/null; then
                    efetch -db nucleotide -id "${ACC}" -format fasta \
                        > "${OUT_FASTA}"
                else
                    # Last resort: plain curl against NCBI eutils REST API
                    log "  (using curl fallback for ${ACC})"
                    curl -fsSL \
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=fasta&retmode=text" \
                        > "${OUT_FASTA}"
                fi
            fi

            [ -s "${OUT_FASTA}" ] || die "Download failed for ${ACC}"
            log "  ${ACC} — saved to ${OUT_FASTA}"
        fi

        FETCHED_FILES+=("${OUT_FASTA}")
    done

    # Concatenate all chromosomes into one multi-chromosome FASTA
    log "Merging ${#FETCHED_FILES[@]} chromosome(s) into ${MERGED_REF} …"
    cat "${FETCHED_FILES[@]}" > "${MERGED_REF}"
    ;;

# ── LOCAL: use an existing merged FASTA as-is ──────────────────────────────
local)
    log "Using local reference: ${REF_PATH}"
    [ -f "${REF_PATH}" ] || die "Reference file not found: ${REF_PATH}"
    cp -n "${REF_PATH}" "${MERGED_REF}" 2>/dev/null || true
    MERGED_REF="${REF_PATH}"
    ;;

# ── CAT: concatenate user-supplied per-chromosome FASTAs ───────────────────
cat)
    log "Concatenating chromosome FASTAs: ${REF_FILES}"
    # Convert space-separated string to an array so paths with spaces work
    # when the user passes REF_FILES as a bash array assignment instead.
    read -ra _REF_ARRAY <<< "${REF_FILES}"
    for f in "${_REF_ARRAY[@]}"; do
        [ -f "$f" ] || die "Reference file not found: $f"
    done
    cat "${_REF_ARRAY[@]}" > "${MERGED_REF}"
    log "Merged reference written to ${MERGED_REF}"
    ;;

*)
    die "Unknown REF_SOURCE '${REF_SOURCE}'. Choose: ncbi | local | cat"
    ;;
esac

# Report how many sequences ended up in the reference
N_CHROMS=$(grep -c "^>" "${MERGED_REF}" 2>/dev/null || echo "0")
[ -z "${N_CHROMS}" ] && N_CHROMS=0
log "Reference ready: ${MERGED_REF}  (${N_CHROMS} sequence(s) / chromosome(s))"

if [ "${N_CHROMS}" -gt 1 ]; then
    log "Multi-chromosome mitogenome detected (${N_CHROMS} chromosomes)."
    log "Per-chromosome baselines, NUMT thresholds and coverage budgets will"
    log "be computed independently for each chromosome."
fi

# ---------------------------------------------------------------------------
# ⑤ BUILD PIPELINE COMMAND
# ---------------------------------------------------------------------------

PIPELINE_CMD=(
    python "${SCRIPT_DIR}/pipeline.py"
    --reads        "${READS}"
    --reference    "${MERGED_REF}"
    --outdir       "${OUTDIR}"
    --threads      "${THREADS}"
    --target-coverage "${TARGET_COV}"
    --genome-size  "${GENOME_SIZE}"
)

if [ "${RUN_ASSEMBLER}" = "1" ]; then
    PIPELINE_CMD+=(--run-assembler --assembler flye)
fi

# ---------------------------------------------------------------------------
# ⑥ RUN
# ---------------------------------------------------------------------------

log "========================================"
log "Starting MitoAssembler pipeline"
log "  Reads      : ${READS}"
log "  Reference  : ${MERGED_REF} (${N_CHROMS} chr)"
log "  Output dir : ${OUTDIR}"
log "  Threads    : ${THREADS}"
log "  Target cov : ${TARGET_COV}x"
log "  Assembler  : $([ "${RUN_ASSEMBLER}" = "1" ] && echo "Flye (${GENOME_SIZE})" || echo "disabled")"
log "========================================"

"${PIPELINE_CMD[@]}"

log "========================================"
log "Pipeline finished."
log "HTML report: ${OUTDIR}/report/report.html"
log "========================================"
