#!/usr/bin/env bash
#
# PalaeoKin.sh — aDNA Relatedness Pipeline using ANGSD + NgsRelate
#
# Tests whether ancient DNA samples originate from the same or different individuals
# using genotype-likelihood methods appropriate for ultra-low coverage aDNA data.
#
# Works with any reference assembly (scaffold-level, chromosome-level, or mixed)
# and any number of BAM files.
#
# Usage:
#   ./PalaeoKin.sh --bams bam_list.txt --ref reference.fna --outdir ./results
#
# See --help for all options.
#
set -euo pipefail

# ============================================================================
# Pipeline version
# ============================================================================
PIPELINE_VERSION="1.0.0"

# ============================================================================
# Resolve script locations
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${SCRIPT_DIR}/scripts"
UTILS="${SCRIPTS}/utils"

# ============================================================================
# Default parameters
ENV_NAME="palaeokin"
# ============================================================================
BAMS=""
REF=""
OUTDIR=""
THREADS=4
DOWNSAMPLE=false
MIXTURE_CHECK=false
TRIM=2
SNP_PVAL="1e-6"
MIN_MAF="0.05"
MIN_CONTIG_SIZE=1000000
FORCE=false

# ============================================================================
# Usage / help
# ============================================================================
usage() {
    cat <<EOF
PalaeoKin v${PIPELINE_VERSION}
===========================================
Tests whether aDNA samples are from the same or different individuals using
genotype-likelihood methods (ANGSD + NgsRelate). Designed for ultra-low coverage
ancient DNA, including cross-species mapping scenarios.

Usage:
  $(basename "$0") --bams FILE --ref FILE --outdir DIR [OPTIONS]

Required arguments:
  --bams FILE              Text file with one BAM path per line
  --ref FILE               Reference FASTA (any assembly level)
  --outdir DIR             Output directory

Optional arguments:
  --threads N              Number of threads (default: ${THREADS})
  --downsample             Enable downsampling validation
  --mixture-check          Enable per-sample mixture diagnostics
  --trim N                 Base pairs to trim from read ends (default: ${TRIM})
  --snp-pval VAL           ANGSD SNP p-value threshold (default: ${SNP_PVAL})
  --min-maf VAL            Minimum minor allele frequency (default: ${MIN_MAF})
  --min-contig-size N      Minimum contig size in bp to include (default: ${MIN_CONTIG_SIZE})
  --force                  Re-run all steps (ignore sentinel files)
  --help                   Show this help message

Examples:
  # Basic run with two samples
  $(basename "$0") --bams samples.txt --ref ref.fna --outdir results

  # Full analysis with validation
  $(basename "$0") --bams samples.txt --ref ref.fna --outdir results \\
      --threads 8 --downsample --mixture-check

  # Fragmented assembly (lower contig size threshold)
  $(basename "$0") --bams samples.txt --ref ref.fna --outdir results \\
      --min-contig-size 100000

Notes:
  - BAM files must be coordinate-sorted and deduplicated
  - Transversion-only SNPs are used to avoid aDNA deamination artifacts
  - Cross-species mapping is supported (uses -C 50 in ANGSD)
  - Pipeline creates/uses a conda environment 'palaeokin' on first run
EOF
    exit 0
}

# ============================================================================
# Logging helpers
# ============================================================================
log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

# ============================================================================
# Parse command-line arguments
# ============================================================================
if [[ $# -eq 0 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bams)        BAMS="$2"; shift 2 ;;
        --ref)         REF="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --downsample)  DOWNSAMPLE=true; shift ;;
        --mixture-check) MIXTURE_CHECK=true; shift ;;
        --trim)        TRIM="$2"; shift 2 ;;
        --snp-pval)    SNP_PVAL="$2"; shift 2 ;;
        --min-maf)     MIN_MAF="$2"; shift 2 ;;
        --min-contig-size) MIN_CONTIG_SIZE="$2"; shift 2 ;;
        --force)       FORCE=true; shift ;;
        --help|-h)     usage ;;
        *)             log_error "Unknown option: $1"; echo "Use --help for usage."; exit 1 ;;
    esac
done

# ============================================================================
# Validate required arguments
# ============================================================================
errors=0
if [[ -z "$BAMS" ]]; then
    log_error "Missing required argument: --bams"
    errors=1
fi
if [[ -z "$REF" ]]; then
    log_error "Missing required argument: --ref"
    errors=1
fi
if [[ -z "$OUTDIR" ]]; then
    log_error "Missing required argument: --outdir"
    errors=1
fi
if [[ $errors -ne 0 ]]; then
    echo "Use --help for usage."
    exit 1
fi

# Resolve to absolute paths
BAMS="$(cd "$(dirname "$BAMS")" && pwd)/$(basename "$BAMS")"
REF="$(cd "$(dirname "$REF")" && pwd)/$(basename "$REF")"

if [[ ! -f "$BAMS" ]]; then
    log_error "BAM list file not found: $BAMS"
    exit 1
fi
if [[ ! -f "$REF" ]]; then
    log_error "Reference FASTA not found: $REF"
    exit 1
fi

# Count samples
N_SAMPLES=$(grep -c . "$BAMS" || true)
if [[ "$N_SAMPLES" -lt 2 ]]; then
    log_error "At least 2 BAM files required. Found ${N_SAMPLES} in ${BAMS}."
    exit 1
fi

# ============================================================================
# Create output directory structure
# ============================================================================
mkdir -p "${OUTDIR}"/{logs,angsd,ngsrelate,validation,reports,plots}
OUTDIR="$(cd "$OUTDIR" && pwd)"
LOGDIR="${OUTDIR}/logs"

# ============================================================================
# Export all parameters for subscripts
# ============================================================================
export PIPELINE_VERSION SCRIPT_DIR SCRIPTS UTILS
export BAMS REF OUTDIR THREADS DOWNSAMPLE MIXTURE_CHECK
export TRIM SNP_PVAL MIN_MAF MIN_CONTIG_SIZE FORCE
export N_SAMPLES LOGDIR

# ============================================================================
# Sentinel file helpers
# ============================================================================
step_done() {
    local step="$1"
    if [[ "$FORCE" == "true" ]]; then
        return 1
    fi
    [[ -f "${OUTDIR}/.${step}.done" ]]
}

mark_done() {
    local step="$1"
    date "+%Y-%m-%d %H:%M:%S" > "${OUTDIR}/.${step}.done"
}

# Export helpers for subscripts
export -f step_done mark_done

# ============================================================================
# Pipeline banner
# ============================================================================
log_info "============================================================"
log_info "  PalaeoKin v${PIPELINE_VERSION}"
log_info "============================================================"
log_info "BAM list:          ${BAMS} (${N_SAMPLES} samples)"
log_info "Reference:         ${REF}"
log_info "Output directory:  ${OUTDIR}"
log_info "Threads:           ${THREADS}"
log_info "Trim:              ${TRIM} bp"
log_info "SNP p-value:       ${SNP_PVAL}"
log_info "Min MAF:           ${MIN_MAF}"
log_info "Min contig size:   ${MIN_CONTIG_SIZE} bp"
log_info "Downsampling:      ${DOWNSAMPLE}"
log_info "Mixture check:     ${MIXTURE_CHECK}"
log_info "Force re-run:      ${FORCE}"
log_info "============================================================"

# ============================================================================
# Step 0: Environment Setup
# ============================================================================
log_info "Step 0: Environment setup"
if step_done "00_setup"; then
    log_info "  -> Skipping (already done). Use --force to re-run."
else
    bash "${SCRIPTS}/00_setup_environment.sh" 2>&1 | tee "${LOGDIR}/00_setup_environment.log"
    mark_done "00_setup"
fi

# Activate conda environment
CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [[ -z "$CONDA_BASE" ]]; then
    log_error "conda not found. Please install conda/miniforge3."
    exit 1
fi
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
log_info "  -> Conda env '${ENV_NAME}' activated"

# ============================================================================
# Step 1: Input Validation
# ============================================================================
log_info "Step 1: Input validation"
if step_done "01_validate"; then
    log_info "  -> Skipping (already done). Use --force to re-run."
else
    bash "${SCRIPTS}/01_validate_inputs.sh" 2>&1 | tee "${LOGDIR}/01_validate_inputs.log"
    mark_done "01_validate"
fi

# ============================================================================
# Step 2: ANGSD Genotype Likelihoods
# ============================================================================
log_info "Step 2: ANGSD genotype likelihoods"
if step_done "02_angsd"; then
    log_info "  -> Skipping (already done). Use --force to re-run."
else
    bash "${SCRIPTS}/02_angsd_genotype_likelihoods.sh" 2>&1 | tee "${LOGDIR}/02_angsd_gl.log"
    mark_done "02_angsd"
fi

# ============================================================================
# Step 3: NgsRelate
# ============================================================================
log_info "Step 3: NgsRelate"
if step_done "03_ngsrelate"; then
    log_info "  -> Skipping (already done). Use --force to re-run."
else
    bash "${SCRIPTS}/03_run_ngsrelate.sh" 2>&1 | tee "${LOGDIR}/03_ngsrelate.log"
    mark_done "03_ngsrelate"
fi

# ============================================================================
# Step 4: Downsampling Validation (optional)
# ============================================================================
if [[ "$DOWNSAMPLE" == "true" ]]; then
    log_info "Step 4: Downsampling validation"
    if step_done "04_downsample"; then
        log_info "  -> Skipping (already done). Use --force to re-run."
    else
        bash "${SCRIPTS}/04_downsample_validation.sh" 2>&1 | tee "${LOGDIR}/04_downsample.log"
        mark_done "04_downsample"
    fi
else
    log_info "Step 4: Downsampling validation (skipped — use --downsample to enable)"
fi

# ============================================================================
# Step 5: Mixture Diagnostics (optional)
# ============================================================================
if [[ "$MIXTURE_CHECK" == "true" ]]; then
    log_info "Step 5: Mixture diagnostics"
    if step_done "05_mixture"; then
        log_info "  -> Skipping (already done). Use --force to re-run."
    else
        bash "${SCRIPTS}/05_mixture_diagnostics.sh" 2>&1 | tee "${LOGDIR}/05_mixture.log"
        mark_done "05_mixture"
    fi
else
    log_info "Step 5: Mixture diagnostics (skipped — use --mixture-check to enable)"
fi

# ============================================================================
# Step 6: Summary Report
# ============================================================================
log_info "Step 6: Generating summary report"
python3 "${UTILS}/summarize_results.py" \
    --outdir "${OUTDIR}" \
    --bams "${BAMS}" \
    --ref "${REF}" \
    --n-samples "${N_SAMPLES}" \
    --downsample "${DOWNSAMPLE}" \
    --mixture-check "${MIXTURE_CHECK}" \
    2>&1 | tee "${LOGDIR}/06_summary.log"

log_info "============================================================"
log_info "  Pipeline complete!"
log_info "============================================================"
log_info "Summary report: ${OUTDIR}/reports/summary_report.txt"
log_info "Plots:          ${OUTDIR}/plots/"
log_info "============================================================"
