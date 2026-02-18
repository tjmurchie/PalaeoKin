#!/usr/bin/env bash
#
# 03_run_ngsrelate.sh — Run NgsRelate and parse results
#
# Extracts allele frequencies from ANGSD output, runs NgsRelate,
# and invokes the Python parser to classify relationships.
#
set -euo pipefail

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

ANGSD_DIR="${OUTDIR}/angsd"
NGSRELATE_DIR="${OUTDIR}/ngsrelate"
PREFIX="${ANGSD_DIR}/joint_gl"

# ============================================================================
# Get NgsRelate binary name
# ============================================================================
if [[ -f "${OUTDIR}/.ngsrelate_bin" ]]; then
    NGSRELATE_BIN=$(cat "${OUTDIR}/.ngsrelate_bin")
else
    # Auto-detect
    if command -v ngsRelate &>/dev/null; then
        NGSRELATE_BIN="ngsRelate"
    elif command -v NgsRelate &>/dev/null; then
        NGSRELATE_BIN="NgsRelate"
    elif command -v ngsrelate &>/dev/null; then
        NGSRELATE_BIN="ngsrelate"
    else
        log_error "NgsRelate binary not found"
        exit 1
    fi
fi
log_info "Using NgsRelate binary: ${NGSRELATE_BIN}"

# ============================================================================
# Extract allele frequencies from ANGSD .mafs.gz
# ============================================================================
log_info "Extracting allele frequencies from ANGSD output..."

# Dynamically find the frequency column (contains "knownEM" or just column with EM frequencies)
# ANGSD .mafs.gz header: chromo position major minor ref knownEM nInd
FREQ_COL=$(zcat "${PREFIX}.mafs.gz" | head -1 | tr '\t' '\n' | grep -n -i "EM" | head -1 | cut -d: -f1)

if [[ -z "$FREQ_COL" ]]; then
    log_warn "Could not find EM frequency column by name. Using column 6 (default ANGSD layout)."
    FREQ_COL=6
fi
log_info "  Frequency column index: ${FREQ_COL}"

# Extract frequencies (skip header)
FREQ_FILE="${NGSRELATE_DIR}/allele_freqs.txt"
zcat "${PREFIX}.mafs.gz" | tail -n +2 | cut -f"${FREQ_COL}" > "${FREQ_FILE}"

N_FREQ=$(wc -l < "${FREQ_FILE}")
log_info "  Extracted ${N_FREQ} frequency values"

# ============================================================================
# Extract sample names from BAM @RG headers or filenames
# ============================================================================
log_info "Extracting sample names..."

NAMES_FILE="${NGSRELATE_DIR}/sample_names.txt"
WORKING_BAMS="${OUTDIR}/working_bam_list.txt"

> "${NAMES_FILE}"
idx=0
while IFS= read -r bam; do
    [[ -z "$bam" ]] && continue
    # Try @RG SM tag first
    rg_name=$(samtools view -H "$bam" 2>/dev/null | grep "^@RG" | head -1 | \
              sed -n 's/.*SM:\([^\t]*\).*/\1/p' || true)
    if [[ -n "$rg_name" ]]; then
        echo "${rg_name}" >> "${NAMES_FILE}"
    else
        # Fall back to filename
        echo "$(basename "$bam" .bam)" >> "${NAMES_FILE}"
    fi
    idx=$((idx + 1))
done < "${WORKING_BAMS}"

log_info "  Sample names:"
while IFS= read -r name; do
    log_info "    - ${name}"
done < "${NAMES_FILE}"

# ============================================================================
# Run NgsRelate
# ============================================================================
log_info "Running NgsRelate..."

N_SITES=$(cat "${PREFIX}.nsites" 2>/dev/null || echo "0")
NGSRELATE_OUT="${NGSRELATE_DIR}/ngsrelate_output.tsv"

"${NGSRELATE_BIN}" \
    -g "${PREFIX}.glf.gz" \
    -f "${FREQ_FILE}" \
    -n "${N_SAMPLES}" \
    -z "${NAMES_FILE}" \
    -O "${NGSRELATE_OUT}" \
    2>&1

if [[ ! -f "${NGSRELATE_OUT}" ]] || [[ ! -s "${NGSRELATE_OUT}" ]]; then
    log_error "NgsRelate produced no output"
    exit 1
fi

log_info "NgsRelate output: ${NGSRELATE_OUT}"

# Show raw output
log_info "NgsRelate raw results:"
cat "${NGSRELATE_OUT}" | while IFS= read -r line; do
    log_info "  ${line}"
done

# ============================================================================
# Parse results with Python classifier
# ============================================================================
log_info "Parsing NgsRelate results..."

python3 "${UTILS}/parse_ngsrelate.py" \
    --input "${NGSRELATE_OUT}" \
    --n-sites "${N_SITES}" \
    --param-set "$(cat "${ANGSD_DIR}/param_set_used.txt" 2>/dev/null || echo "standard")" \
    --outdir "${NGSRELATE_DIR}"

log_info "NgsRelate analysis complete."
log_info "Results: ${NGSRELATE_DIR}/relatedness_results.tsv"
log_info "Report:  ${NGSRELATE_DIR}/relatedness_report.txt"
