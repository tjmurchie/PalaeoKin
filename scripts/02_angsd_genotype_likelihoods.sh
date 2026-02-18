#!/usr/bin/env bash
#
# 02_angsd_genotype_likelihoods.sh — ANGSD joint genotype likelihoods at TV sites
#
# Computes genotype likelihoods for all samples jointly using ANGSD,
# restricted to transversion SNPs to avoid aDNA deamination artifacts.
# Includes adaptive parameter relaxation if initial run yields too few sites.
#
set -euo pipefail

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

ANGSD_DIR="${OUTDIR}/angsd"
WORKING_BAMS="${OUTDIR}/working_bam_list.txt"
REGIONS_FILE="${ANGSD_DIR}/regions.txt"
PREFIX="${ANGSD_DIR}/joint_gl"

# ============================================================================
# Run ANGSD with standard parameters
# ============================================================================
run_angsd() {
    local prefix="$1"
    local snp_pval="$2"
    local min_maf="$3"
    local trim="$4"
    local label="$5"

    log_info "Running ANGSD (${label})..."
    log_info "  SNP_pval=${snp_pval}, minMaf=${min_maf}, trim=${trim}"

    angsd \
        -bam "${WORKING_BAMS}" \
        -ref "${REF}" \
        -out "${prefix}" \
        -GL 1 \
        -doGlf 3 \
        -doMajorMinor 1 \
        -doMaf 1 \
        -SNP_pval "${snp_pval}" \
        -minMaf "${min_maf}" \
        -rmTrans 1 \
        -minMapQ 30 \
        -minQ 20 \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 0 \
        -C 50 \
        -doCounts 1 \
        -trim "${trim}" \
        -minInd "${N_SAMPLES}" \
        -setMinDepth 2 \
        -rf "${REGIONS_FILE}" \
        -nThreads "${THREADS}" \
        2>&1

    # Count sites
    if [[ -f "${prefix}.mafs.gz" ]]; then
        local n_sites
        n_sites=$(zcat "${prefix}.mafs.gz" | tail -n +2 | wc -l)
        echo "${n_sites}" > "${prefix}.nsites"
        log_info "  Sites found: ${n_sites}"
        return 0
    else
        log_warn "  ANGSD produced no output"
        echo "0" > "${prefix}.nsites"
        return 1
    fi
}

# ============================================================================
# Verify transversion-only output
# ============================================================================
verify_transversions() {
    local prefix="$1"

    log_info "Verifying transversion-only SNPs..."

    # Check that no transitions leaked through
    # Transitions: A<>G, C<>T
    local n_transitions
    n_transitions=$(zcat "${prefix}.mafs.gz" | tail -n +2 | awk '
        ($3 == "A" && $4 == "G") || ($3 == "G" && $4 == "A") ||
        ($3 == "C" && $4 == "T") || ($3 == "T" && $4 == "C") {count++}
        END {print count+0}
    ')

    if [[ "$n_transitions" -gt 0 ]]; then
        log_warn "  Found ${n_transitions} transitions in output — -rmTrans may not have worked."
        log_warn "  Post-filtering to remove transitions..."

        # Create a sites file with only transversions
        local tv_sites="${ANGSD_DIR}/tv_sites.txt"
        zcat "${prefix}.mafs.gz" | tail -n +2 | awk '
            !(($3 == "A" && $4 == "G") || ($3 == "G" && $4 == "A") ||
              ($3 == "C" && $4 == "T") || ($3 == "T" && $4 == "C")) {
                print $1"\t"$2
            }
        ' > "${tv_sites}"

        local n_tv
        n_tv=$(wc -l < "${tv_sites}")
        log_info "  Transversion sites after filtering: ${n_tv}"

        if [[ "$n_tv" -gt 0 ]]; then
            # Index the sites file for ANGSD
            angsd sites index "${tv_sites}"

            # Re-run ANGSD with sites file
            log_info "  Re-running ANGSD with transversion-only sites file..."
            angsd \
                -bam "${WORKING_BAMS}" \
                -ref "${REF}" \
                -out "${prefix}" \
                -GL 1 \
                -doGlf 3 \
                -doMajorMinor 1 \
                -doMaf 1 \
                -minMapQ 30 \
                -minQ 20 \
                -uniqueOnly 1 \
                -remove_bads 1 \
                -only_proper_pairs 0 \
                -C 50 \
                -doCounts 1 \
                -trim "${TRIM}" \
                -minInd "${N_SAMPLES}" \
                -setMinDepth 2 \
                -sites "${tv_sites}" \
                -rf "${REGIONS_FILE}" \
                -nThreads "${THREADS}" \
                2>&1

            n_sites=$(zcat "${prefix}.mafs.gz" | tail -n +2 | wc -l)
            echo "${n_sites}" > "${prefix}.nsites"
            log_info "  Final site count after TV filtering: ${n_sites}"
        fi
    else
        log_info "  Confirmed: all sites are transversions."
    fi
}

# ============================================================================
# Main execution
# ============================================================================

# Run with standard parameters
run_angsd "${PREFIX}" "${SNP_PVAL}" "${MIN_MAF}" "${TRIM}" "standard parameters"

N_SITES=$(cat "${PREFIX}.nsites" 2>/dev/null || echo "0")

# Track which parameter set was used
PARAM_SET="standard"

# ============================================================================
# Adaptive relaxation if too few sites
# ============================================================================
if [[ "$N_SITES" -lt 50 ]]; then
    log_warn "Only ${N_SITES} sites found — attempting adaptive relaxation..."

    # Level 1: relax p-value and MAF
    RELAXED_PREFIX="${ANGSD_DIR}/joint_gl_relaxed1"
    run_angsd "${RELAXED_PREFIX}" "1e-4" "0.01" "${TRIM}" "relaxed p-val + MAF"
    N_RELAXED1=$(cat "${RELAXED_PREFIX}.nsites" 2>/dev/null || echo "0")

    if [[ "$N_RELAXED1" -gt "$N_SITES" ]]; then
        log_info "Relaxed parameters yielded ${N_RELAXED1} sites (was ${N_SITES})"
        cp "${RELAXED_PREFIX}.mafs.gz" "${PREFIX}.mafs.gz"
        cp "${RELAXED_PREFIX}.glf.gz" "${PREFIX}.glf.gz"
        echo "${N_RELAXED1}" > "${PREFIX}.nsites"
        N_SITES=$N_RELAXED1
        PARAM_SET="relaxed_pval_maf"
    fi

    # Level 2: also remove trimming
    if [[ "$N_SITES" -lt 50 ]]; then
        RELAXED_PREFIX2="${ANGSD_DIR}/joint_gl_relaxed2"
        run_angsd "${RELAXED_PREFIX2}" "1e-4" "0.01" "0" "relaxed p-val + MAF + no trim"
        N_RELAXED2=$(cat "${RELAXED_PREFIX2}.nsites" 2>/dev/null || echo "0")

        if [[ "$N_RELAXED2" -gt "$N_SITES" ]]; then
            log_info "Further relaxation yielded ${N_RELAXED2} sites (was ${N_SITES})"
            cp "${RELAXED_PREFIX2}.mafs.gz" "${PREFIX}.mafs.gz"
            cp "${RELAXED_PREFIX2}.glf.gz" "${PREFIX}.glf.gz"
            echo "${N_RELAXED2}" > "${PREFIX}.nsites"
            N_SITES=$N_RELAXED2
            PARAM_SET="relaxed_all"
        fi
    fi
fi

# Save which parameter set was actually used
echo "${PARAM_SET}" > "${ANGSD_DIR}/param_set_used.txt"

if [[ "$PARAM_SET" != "standard" ]]; then
    log_warn "*** Results use RELAXED parameters (${PARAM_SET}) — interpret with extra caution ***"
fi

# ============================================================================
# Verify transversions
# ============================================================================
if [[ "$N_SITES" -gt 0 ]]; then
    verify_transversions "${PREFIX}"
    # Update site count after verification
    N_SITES=$(zcat "${PREFIX}.mafs.gz" | tail -n +2 | wc -l)
    echo "${N_SITES}" > "${PREFIX}.nsites"
fi

# ============================================================================
# Final status
# ============================================================================
if [[ "$N_SITES" -eq 0 ]]; then
    log_error "No informative transversion sites found. Cannot proceed with relatedness estimation."
    log_error "Possible causes:"
    log_error "  - Coverage too low for overlap between samples"
    log_error "  - Reference too divergent from samples"
    log_error "  - BAM files have no mapped reads to large contigs"
    exit 1
fi

log_info "ANGSD complete: ${N_SITES} transversion sites (parameters: ${PARAM_SET})"
log_info "Output files:"
log_info "  GLs:  ${PREFIX}.glf.gz"
log_info "  MAFs: ${PREFIX}.mafs.gz"
