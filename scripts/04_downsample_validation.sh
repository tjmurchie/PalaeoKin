#!/usr/bin/env bash
#
# 04_downsample_validation.sh — Robustness validation via downsampling
#
# Downsamples the best BAM to several fractions (50%, 25%, 10%, and
# matched-to-worst) and re-runs ANGSD + NgsRelate at each level.
# Compares KING estimates across levels for consistency.
#
set -euo pipefail

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

VALIDATION_DIR="${OUTDIR}/validation/downsampling"
ANGSD_DIR="${OUTDIR}/angsd"
REGIONS_FILE="${ANGSD_DIR}/regions.txt"
WORKING_BAMS="${OUTDIR}/working_bam_list.txt"

mkdir -p "${VALIDATION_DIR}"

# ----------------------------------------------------------------------------
# Helper: robustly extract KING/R0/R1/nSites from an NgsRelate TSV (v2)
# ----------------------------------------------------------------------------
extract_ngsrelate_stats() {
    local tsv="$1"
    if [[ ! -s "${tsv}" ]]; then
        echo -e "NA\tNA\tNA\tNA"
        return 0
    fi

    awk -F'\t' '
        NR==1{
            for(i=1;i<=NF;i++){
                if($i=="KING") king=i
                else if($i=="R0") r0=i
                else if($i=="R1") r1=i
                else if($i=="nSites") ns=i
            }
            next
        }
        NR==2{
            printf "%s\t%s\t%s\t%s\n", (king? $king:"NA"), (r0? $r0:"NA"), (r1? $r1:"NA"), (ns? $ns:"NA")
            exit
        }' "${tsv}"
}

# ============================================================================
# Identify best and worst samples
# ============================================================================
BEST_SAMPLE=$(cat "${OUTDIR}/validation/best_sample.txt")
WORST_SAMPLE=$(cat "${OUTDIR}/validation/worst_sample.txt")
BEST_READS=$(cat "${OUTDIR}/validation/best_reads.txt")
WORST_READS=$(cat "${OUTDIR}/validation/worst_reads.txt")

# Find the actual BAM paths
BEST_BAM=$(grep "${BEST_SAMPLE}" "${WORKING_BAMS}" | head -1)
WORST_BAM=$(grep "${WORST_SAMPLE}" "${WORKING_BAMS}" | head -1)

if [[ -z "$BEST_BAM" ]] || [[ -z "$WORST_BAM" ]]; then
    log_error "Could not identify best/worst BAM files"
    exit 1
fi

log_info "Best BAM:  ${BEST_BAM} (${BEST_READS} reads)"
log_info "Worst BAM: ${WORST_BAM} (${WORST_READS} reads)"

# ============================================================================
# Calculate downsampling fractions
# ============================================================================
# Match worst: downsample best to same read count as worst
MATCH_FRAC=$(awk "BEGIN {f = ${WORST_READS} / ${BEST_READS}; if(f > 1) f = 1; printf \"%.6f\", f}")
log_info "Match-worst fraction: ${MATCH_FRAC}"

FRACTIONS=("0.50" "0.25" "0.10" "${MATCH_FRAC}")
LABELS=("50pct" "25pct" "10pct" "matched")

# ============================================================================
# NgsRelate binary (auto-detect, same logic as 03_run_ngsrelate.sh)
# ============================================================================
if [[ -f "${OUTDIR}/.ngsrelate_bin" ]]; then
    NGSRELATE_BIN=$(cat "${OUTDIR}/.ngsrelate_bin")
else
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
# Results collector
# ============================================================================
RESULTS_FILE="${VALIDATION_DIR}/downsampling_results.tsv"
echo -e "fraction\tlabel\tn_sites\tking\tR0\tR1\tnSNPs_ngsrelate" > "${RESULTS_FILE}"

# Add full-data result (parsed robustly by header)
FULL_NGSRELATE="${OUTDIR}/ngsrelate/ngsrelate_output.tsv"
if [[ -s "$FULL_NGSRELATE" ]]; then
    FULL_NSITES=$(cat "${ANGSD_DIR}/joint_gl.nsites" 2>/dev/null || echo "0")
    read -r FULL_KING FULL_R0 FULL_R1 FULL_NSNP <<<"$(extract_ngsrelate_stats "${FULL_NGSRELATE}")"
    echo -e "1.00\tfull\t${FULL_NSITES}\t${FULL_KING}\t${FULL_R0}\t${FULL_R1}\t${FULL_NSNP}" >> "${RESULTS_FILE}"
else
    log_warn "Full NgsRelate output not found at ${FULL_NGSRELATE}; downsampling will still run, but no full baseline row will be added."
fi

# ============================================================================
# Run downsampling iterations
# ============================================================================
for i in "${!FRACTIONS[@]}"; do
    frac="${FRACTIONS[$i]}"
    label="${LABELS[$i]}"

    log_info "=========================================="
    log_info "Downsampling: ${label} (fraction=${frac})"
    log_info "=========================================="

    DS_DIR="${VALIDATION_DIR}/${label}"
    mkdir -p "${DS_DIR}"

    # Downsample best BAM
    DS_BAM="${DS_DIR}/$(basename "$BEST_BAM" .bam).ds_${label}.bam"

    # Handle fraction ~1.0 (no downsampling; avoids invalid samtools -s seed.fraction strings)
    IS_ONE=$(awk -v f="${frac}" 'BEGIN{ if(f>=0.999999) print 1; else print 0 }')
    if [[ "${IS_ONE}" -eq 1 ]]; then
        log_info "  Fraction ~1.0; copying BAM without downsampling..."
        samtools view -b "${BEST_BAM}" > "${DS_BAM}"
    else
        # samtools view -s uses SEED.FRACTION format
        SEED=42
        SUBSAMPLE="${SEED}.${frac#0.}"  # e.g., 42.50 for 50%
        log_info "  Downsampling with samtools view -s ${SUBSAMPLE}..."
        samtools view -b -s "${SUBSAMPLE}" "${BEST_BAM}" > "${DS_BAM}"
    fi
    samtools index "${DS_BAM}"

    # Create BAM list for this iteration (downsampled best + all other BAMs unchanged)
    DS_BAMLIST="${DS_DIR}/bam_list.txt"
    echo "${DS_BAM}" > "${DS_BAMLIST}"
    while IFS= read -r bam; do
        [[ -z "$bam" ]] && continue
        if [[ "$(readlink -f "$bam")" != "$(readlink -f "$BEST_BAM")" ]]; then
            echo "$bam" >> "${DS_BAMLIST}"
        fi
    done < "${WORKING_BAMS}"

    DS_N=$(wc -l < "${DS_BAMLIST}")
    DS_PREFIX="${DS_DIR}/angsd_ds"

    # Count reads in downsampled BAM
    DS_READS=$(samtools idxstats "${DS_BAM}" | awk '{sum += $3} END {print sum}')
    log_info "  Downsampled reads: ${DS_READS}"

    # Run ANGSD (same core params as main joint run)
    log_info "  Running ANGSD..."
    set +e
    angsd         -bam "${DS_BAMLIST}"         -ref "${REF}"         -out "${DS_PREFIX}"         -GL 1         -doGlf 3         -doMajorMinor 1         -doMaf 1         -SNP_pval "${SNP_PVAL}"         -minMaf "${MIN_MAF}"         -rmTrans 1         -minMapQ 30         -minQ 20         -uniqueOnly 1         -remove_bads 1         -only_proper_pairs 0         -C 50         -doCounts 1         -trim "${TRIM}"         -minInd "${DS_N}"         -setMinDepth 2         -rf "${REGIONS_FILE}"         -nThreads "${THREADS}"         2>&1
    ANGSD_EXIT=$?
    set -e
    if [[ "${ANGSD_EXIT}" -ne 0 ]]; then
        log_warn "  ANGSD returned non-zero exit code (${ANGSD_EXIT}) for ${label}. Continuing."
    fi

    # Count sites
    DS_SITES=0
    if [[ -f "${DS_PREFIX}.mafs.gz" ]]; then
        DS_SITES=$(zcat "${DS_PREFIX}.mafs.gz" | tail -n +2 | wc -l)
    fi
    log_info "  Sites found: ${DS_SITES}"

    # Run NgsRelate if we have sites
    DS_KING="NA"
    DS_R0="NA"
    DS_R1="NA"
    DS_NSNP="NA"

    if [[ "$DS_SITES" -gt 0 ]]; then
        # Extract frequencies
        FREQ_COL=$(zcat "${DS_PREFIX}.mafs.gz" | head -1 | tr '\t' '\n' | grep -n -i "EM" | head -1 | cut -d: -f1)
        [[ -z "$FREQ_COL" ]] && FREQ_COL=6
        zcat "${DS_PREFIX}.mafs.gz" | tail -n +2 | cut -f"${FREQ_COL}" > "${DS_DIR}/freqs.txt"

        # Sample names (same as main run)
        NAMES_FILE="${OUTDIR}/ngsrelate/sample_names.txt"
        DS_NGSRELATE="${DS_DIR}/ngsrelate_output.tsv"

        set +e
        "${NGSRELATE_BIN}"             -g "${DS_PREFIX}.glf.gz"             -f "${DS_DIR}/freqs.txt"             -n "${DS_N}"             -z "${NAMES_FILE}"             -O "${DS_NGSRELATE}"             2>&1
        NGR_EXIT=$?
        set -e
        if [[ "${NGR_EXIT}" -ne 0 ]]; then
            log_warn "  NgsRelate returned non-zero exit code (${NGR_EXIT}) for ${label}. Continuing."
        fi

        if [[ -s "${DS_NGSRELATE}" ]]; then
            read -r DS_KING DS_R0 DS_R1 DS_NSNP <<<"$(extract_ngsrelate_stats "${DS_NGSRELATE}")"
        else
            log_warn "  NgsRelate output missing/empty for ${label}"
        fi
    else
        log_warn "  No sites at ${label} downsampling — skipping NgsRelate"
    fi

    echo -e "${frac}\t${label}\t${DS_SITES}\t${DS_KING}\t${DS_R0}\t${DS_R1}\t${DS_NSNP}" >> "${RESULTS_FILE}"
    log_info "  KING=${DS_KING}, R0=${DS_R0}, R1=${DS_R1}, nSites=${DS_NSNP}"
done

# ============================================================================
# Summary
# ============================================================================
log_info "Downsampling validation complete."
log_info "Results: ${RESULTS_FILE}"
log_info ""
log_info "Downsampling summary:"
column -t -s $'\t' "${RESULTS_FILE}" | while IFS= read -r line; do
    log_info "  ${line}"
done
