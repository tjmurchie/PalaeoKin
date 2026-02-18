#!/usr/bin/env bash
#
# 05_mixture_diagnostics.sh — Per-sample mixture / allele-balance check
#
# Purpose:
#   A lightweight heuristic to flag potential multi-individual mixtures using
#   allele-balance at SNP sites from ANGSD dumpCounts output.
#
# Notes:
#   - This is NOT a formal contamination estimator.
#   - Requires enough mapped reads to produce sites with depth >= 3.
#   - For samples below the threshold, we still emit a report marking them as
#     INSUFFICIENT_READS so the summary/QC panels can include them.
#
set -euo pipefail

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }

MIXTURE_DIR="${OUTDIR}/validation/mixture"
ANGSD_DIR="${OUTDIR}/angsd"
REGIONS_FILE="${ANGSD_DIR}/regions.txt"
WORKING_BAMS="${OUTDIR}/working_bam_list.txt"
STATS_FILE="${OUTDIR}/validation/sample_stats.tsv"

mkdir -p "${MIXTURE_DIR}"

# Minimum reads to attempt ANGSD dumpCounts. Still produces a placeholder report for smaller samples.
MIN_READS_FOR_MIXTURE=${MIN_READS_FOR_MIXTURE:-500000}

# Output formats for plots (comma-separated): png,pdf,svg
PLOT_FORMATS=${PLOT_FORMATS:-"png,pdf,svg"}

SAMPLES_TOTAL=0
SAMPLES_RUN=0

while IFS= read -r bam; do
    [[ -z "$bam" ]] && continue
    SAMPLES_TOTAL=$((SAMPLES_TOTAL + 1))

    sample_name="$(basename "$bam" .bam)"
    SAMPLE_DIR="${MIXTURE_DIR}/${sample_name}"
    mkdir -p "${SAMPLE_DIR}"
    SAMPLE_PREFIX="${SAMPLE_DIR}/mixture"

    # Get read count from stats (sample_stats.tsv is tab-delimited)
    total_reads=$(awk -v s="${sample_name}" '$1 == s {print $3}' "${STATS_FILE}" | tail -1)
    total_reads=${total_reads:-0}

    if [[ "${total_reads}" -lt "${MIN_READS_FOR_MIXTURE}" ]]; then
        log_info "Mixture check: ${sample_name} has ${total_reads} reads (< ${MIN_READS_FOR_MIXTURE}); writing placeholder report."
        cat > "${SAMPLE_DIR}/mixture_report.txt" <<EOF
Sample: ${sample_name}
Status: INSUFFICIENT_READS
Reads: ${total_reads}
Threshold: ${MIN_READS_FOR_MIXTURE}
Reason: Too few mapped reads to generate enough sites at depth >= 3 for allele-balance diagnostics.
EOF
        continue
    fi

    log_info "Running mixture diagnostics for ${sample_name} (${total_reads} reads)..."

    # Create single-sample BAM list
    echo "$bam" > "${SAMPLE_DIR}/bam_list.txt"

    # Run ANGSD per-sample: get allele counts at transversion sites
    log_info "  Running ANGSD dumpCounts..."
    angsd \
        -bam "${SAMPLE_DIR}/bam_list.txt" \
        -ref "${REF}" \
        -out "${SAMPLE_PREFIX}" \
        -GL 1 \
        -doMajorMinor 1 \
        -doMaf 1 \
        -doCounts 1 \
        -dumpCounts 4 \
        -SNP_pval "1e-4" \
        -rmTrans 1 \
        -minMapQ 30 \
        -minQ 20 \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 0 \
        -C 50 \
        -trim "${TRIM}" \
        -setMinDepth 3 \
        -rf "${REGIONS_FILE}" \
        -nThreads "${THREADS}" \
        2>&1 || true

    if [[ ! -f "${SAMPLE_PREFIX}.mafs.gz" ]] || [[ ! -f "${SAMPLE_PREFIX}.counts.gz" ]]; then
        log_warn "  No ANGSD output for ${sample_name} (likely too few sites with depth >= 3)."
        cat > "${SAMPLE_DIR}/mixture_report.txt" <<EOF
Sample: ${sample_name}
Status: INSUFFICIENT_DATA
Reason: ANGSD did not produce .mafs.gz/.counts.gz (likely too few sites at depth >= 3).
EOF
        continue
    fi

    N_SITES=$(zcat "${SAMPLE_PREFIX}.mafs.gz" | tail -n +2 | wc -l)
    log_info "  Candidate sites (ANGSD mafs): ${N_SITES}"

    if [[ "${N_SITES}" -lt 10 ]]; then
        log_warn "  Too few sites (${N_SITES}) for meaningful mixture analysis."
        cat > "${SAMPLE_DIR}/mixture_report.txt" <<EOF
Sample: ${sample_name}
Status: INSUFFICIENT_DATA
Sites: ${N_SITES}
Reason: Too few SNP sites for allele-balance diagnostics.
EOF
        continue
    fi

    log_info "  Computing allele-balance (counts-based)..."
    python3 "${UTILS}/mixture_analysis.py" \
        --counts "${SAMPLE_PREFIX}.counts.gz" \
        --outdir "${SAMPLE_DIR}" \
        --sample "${sample_name}" \
        --formats "${PLOT_FORMATS}" \
        2>&1

    SAMPLES_RUN=$((SAMPLES_RUN + 1))
done < "${WORKING_BAMS}"

# Consolidate results
SUMMARY="${MIXTURE_DIR}/mixture_summary.txt"
{
  echo "Mixture Diagnostics Summary"
  echo "=========================="
  echo ""
  echo "Notes:"
  echo "  - This is an allele-balance / mixture heuristic (not a formal contamination estimate)."
  echo "  - Low depth causes discrete allele fractions; interpret cautiously."
  echo ""
  echo "Thresholds:"
  echo "  MIN_READS_FOR_MIXTURE = ${MIN_READS_FOR_MIXTURE}"
  echo ""
  echo "Samples processed: ${SAMPLES_TOTAL}"
  echo "Samples with full mixture analysis: ${SAMPLES_RUN}"
  echo ""
  for report in "${MIXTURE_DIR}"/*/mixture_report.txt; do
      [[ -f "$report" ]] || continue
      cat "$report"
      echo ""
      echo "---"
      echo ""
  done
} > "${SUMMARY}"

log_info "Mixture diagnostics complete. Summary: ${SUMMARY}"
