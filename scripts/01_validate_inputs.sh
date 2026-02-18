#!/usr/bin/env bash
#
# 01_validate_inputs.sh — Validate inputs, index files, compute stats
#
# - Indexes reference (.fai) if missing
# - Auto-detects usable contigs by SIZE (not name) — works with any assembly
# - Indexes BAMs if .bai missing (with fallback to outdir symlinks)
# - Computes per-sample stats (reads, avg length, coverage)
# - Estimates expected overlapping TV SNPs and warns if too few
#
set -euo pipefail

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

# ============================================================================
# Index reference genome if needed
# ============================================================================
log_info "Checking reference index..."
if [[ ! -f "${REF}.fai" ]]; then
    log_info "  Creating reference index: ${REF}.fai"
    samtools faidx "$REF"
else
    log_info "  Reference index found: ${REF}.fai"
fi

# ============================================================================
# Auto-detect usable contigs by SIZE
# ============================================================================
log_info "Filtering contigs >= ${MIN_CONTIG_SIZE} bp..."
REGIONS_FILE="${OUTDIR}/angsd/regions.txt"

awk -v min_size="${MIN_CONTIG_SIZE}" '$2 >= min_size {print $1}' "${REF}.fai" \
    > "${REGIONS_FILE}"

N_CONTIGS=$(wc -l < "${REGIONS_FILE}")
TOTAL_BP=$(awk -v min_size="${MIN_CONTIG_SIZE}" '$2 >= min_size {sum += $2} END {print sum+0}' "${REF}.fai")

if [[ "$N_CONTIGS" -eq 0 ]]; then
    log_error "No contigs >= ${MIN_CONTIG_SIZE} bp found in reference."
    log_error "Try lowering --min-contig-size (e.g., --min-contig-size 100000)"
    exit 1
fi

log_info "  Found ${N_CONTIGS} contigs totaling $(numfmt --to=iec ${TOTAL_BP} 2>/dev/null || echo "${TOTAL_BP}") bp"

# ============================================================================
# Index BAM files if needed
# ============================================================================
log_info "Checking BAM files and indexes..."

# We may need to create a working BAM list if we symlink
WORKING_BAMS="${OUTDIR}/working_bam_list.txt"
> "${WORKING_BAMS}"

while IFS= read -r bam; do
    # Skip empty lines
    [[ -z "$bam" ]] && continue

    if [[ ! -f "$bam" ]]; then
        log_error "BAM file not found: $bam"
        exit 1
    fi

    # Check if .bai exists (try both naming conventions)
    bai_found=false
    if [[ -f "${bam}.bai" ]] || [[ -f "${bam%.bam}.bai" ]]; then
        bai_found=true
    fi

    if [[ "$bai_found" == "false" ]]; then
        log_info "  Indexing: $(basename "$bam")"
        # Try indexing in place first
        if samtools index "$bam" 2>/dev/null; then
            bai_found=true
        else
            # If no write permission, symlink to outdir
            log_warn "  Cannot write index to source dir. Symlinking to outdir."
            local_bam="${OUTDIR}/angsd/$(basename "$bam")"
            if [[ ! -L "$local_bam" ]] && [[ ! -f "$local_bam" ]]; then
                ln -s "$(readlink -f "$bam")" "$local_bam"
            fi
            samtools index "$local_bam"
            echo "$local_bam" >> "${WORKING_BAMS}"
            continue
        fi
    fi

    echo "$bam" >> "${WORKING_BAMS}"
done < "$BAMS"

# Verify working BAM list has same count
N_WORKING=$(wc -l < "${WORKING_BAMS}")
if [[ "$N_WORKING" -ne "$N_SAMPLES" ]]; then
    log_error "Working BAM list has ${N_WORKING} entries, expected ${N_SAMPLES}"
    exit 1
fi

# ============================================================================
# Compute per-sample statistics
# ============================================================================
log_info "Computing per-sample statistics..."

STATS_FILE="${OUTDIR}/validation/sample_stats.tsv"
echo -e "sample\tbam_path\ttotal_reads\tavg_read_length\test_coverage\tread_source" > "${STATS_FILE}"

BEST_READS=0
BEST_SAMPLE=""
WORST_READS=999999999999
WORST_SAMPLE=""

while IFS= read -r bam; do
    [[ -z "$bam" ]] && continue
    sample_name="$(basename "$bam" .bam)"

    # Get stats via samtools idxstats (fast, uses index)
    total_reads=$(samtools idxstats "$bam" | awk '{sum += $3} END {print sum}')

    # Get average read length from first 10000 reads via samtools view
    # Note: set +o pipefail in subshell to avoid SIGPIPE exit (head closes pipe early)
    avg_len=$(set +o pipefail; samtools view "$bam" | head -10000 | awk '{sum += length($10); n++} END {if(n>0) printf "%.0f", sum/n; else print 0}')

    # Estimate coverage: (reads * avg_len) / genome_size
    if [[ "$avg_len" -gt 0 ]] && [[ "$TOTAL_BP" -gt 0 ]]; then
        est_cov=$(awk "BEGIN {printf \"%.6f\", (${total_reads} * ${avg_len}) / ${TOTAL_BP}}")
    else
        est_cov="0"
    fi

    echo -e "${sample_name}\t${bam}\t${total_reads}\t${avg_len}\t${est_cov}\t$(basename "$bam")" >> "${STATS_FILE}"

    log_info "  ${sample_name}: ${total_reads} reads, avg ${avg_len}bp, ~${est_cov}x coverage"

    # Track best/worst
    if [[ "$total_reads" -gt "$BEST_READS" ]]; then
        BEST_READS=$total_reads
        BEST_SAMPLE="$sample_name"
    fi
    if [[ "$total_reads" -lt "$WORST_READS" ]]; then
        WORST_READS=$total_reads
        WORST_SAMPLE="$sample_name"
    fi
done < "${WORKING_BAMS}"

log_info "  Best sample:  ${BEST_SAMPLE} (${BEST_READS} reads)"
log_info "  Worst sample: ${WORST_SAMPLE} (${WORST_READS} reads)"

# Save best/worst info for downstream scripts
echo "${BEST_SAMPLE}" > "${OUTDIR}/validation/best_sample.txt"
echo "${WORST_SAMPLE}" > "${OUTDIR}/validation/worst_sample.txt"
echo "${BEST_READS}" > "${OUTDIR}/validation/best_reads.txt"
echo "${WORST_READS}" > "${OUTDIR}/validation/worst_reads.txt"

# ============================================================================
# Coverage feasibility estimate
# ============================================================================
log_info "Estimating expected informative sites..."

# For two samples, the expected number of overlapping TV sites is approximately:
# genome_size * TV_rate * P(covered_sample1) * P(covered_sample2) * P(polymorphic)
# TV_rate ~ 1/3 of all SNPs, divergence ~ 0.001-0.01 for cross-species
# P(covered) ~ 1 - exp(-coverage)

best_cov=$(awk -v s="${BEST_SAMPLE}" '$1 == s {print $5}' "${STATS_FILE}" | tail -1)
worst_cov=$(awk -v s="${WORST_SAMPLE}" '$1 == s {print $5}' "${STATS_FILE}" | tail -1)

# Estimate: total_bp * heterozygosity_rate * tv_fraction * p_cover_best * p_cover_worst
# Using conservative heterozygosity ~ 0.001 and TV fraction ~ 1/3
est_sites=$(python3 -c "
import math
genome = ${TOTAL_BP}
cov_best = ${best_cov}
cov_worst = ${worst_cov}
het_rate = 0.001      # conservative for cross-species
tv_frac = 1.0/3.0     # fraction of SNPs that are transversions
p_best = 1 - math.exp(-cov_best) if cov_best > 0 else 0
p_worst = 1 - math.exp(-cov_worst) if cov_worst > 0 else 0
expected = genome * het_rate * tv_frac * p_best * p_worst
print(int(expected))
" 2>/dev/null || echo "0")

log_info "  Estimated overlapping TV SNPs: ~${est_sites}"

if [[ "$est_sites" -lt 50 ]]; then
    log_warn "  *** VERY LOW expected site count (<50). Results may be unreliable."
    log_warn "  *** Consider: lower --snp-pval, lower --min-maf, or more sequencing."
elif [[ "$est_sites" -lt 100 ]]; then
    log_warn "  *** LOW expected site count (<100). Results should be interpreted with caution."
elif [[ "$est_sites" -lt 500 ]]; then
    log_info "  Marginal but workable site count. Pipeline will use adaptive parameter relaxation if needed."
fi

echo "${est_sites}" > "${OUTDIR}/validation/estimated_sites.txt"

log_info "Input validation complete."
