#!/usr/bin/env bash
#
# 00_setup_environment.sh — Create/verify the PalaeoKin conda environment
#
# Creates the 'palaeokin' conda environment if it doesn't exist, installs required
# dependencies (via environment.yml when available), and verifies they work.
#
set -euo pipefail

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

ENV_NAME="palaeokin"

# Resolve repo root (scripts/..)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_YML="${REPO_ROOT}/environment.yml"

# Check if conda is available
CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [[ -z "${CONDA_BASE}" ]]; then
    log_error "conda/miniforge3 not found. Please install it first."
    exit 1
fi

# Source conda for this shell
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Create environment if it doesn't exist
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    log_info "Conda environment '${ENV_NAME}' already exists."
    log_info "Activating to verify tools..."
    conda activate "${ENV_NAME}"
else
    if [[ -f "${ENV_YML}" ]]; then
        log_info "Creating conda environment '${ENV_NAME}' from environment.yml..."
        # environment.yml must have name: palaeokin
        conda env create -f "${ENV_YML}"
        conda activate "${ENV_NAME}"
    else
        log_warn "environment.yml not found; falling back to explicit conda create."
        conda create -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
            angsd \
            ngsrelate \
            samtools \
            "python>=3.10" \
            numpy \
            pandas \
            matplotlib \
            scipy
        conda activate "${ENV_NAME}"
    fi
fi

# Verify required tools
log_info "Verifying installed tools..."

# ANGSD
if ! command -v angsd &>/dev/null; then
    log_error "angsd not found in environment '${ENV_NAME}'"
    exit 1
fi
angsd -h >/dev/null 2>&1 || true
log_info "  OK: angsd"

# NgsRelate
if ! command -v ngsRelate &>/dev/null; then
    log_error "NgsRelate not found in environment '${ENV_NAME}'"
    exit 1
fi
NgsRelate 2>/dev/null | head -n 1 >/dev/null || true
log_info "  OK: NgsRelate"

# samtools
if ! command -v samtools &>/dev/null; then
    log_error "samtools not found in environment '${ENV_NAME}'"
    exit 1
fi
samtools --version | head -n 1 >/dev/null
log_info "  OK: samtools"

# python deps
python - <<'PY'
import numpy, pandas, matplotlib, scipy
print("OK: python deps")
PY
log_info "  OK: python deps (numpy/pandas/matplotlib/scipy)"

log_info "Environment '${ENV_NAME}' is ready."
