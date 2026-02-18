#!/usr/bin/env python3
"""
summarize_results.py — Generate consolidated summary report and publication-ready plots.

Outputs (all in ./plots/ as PNG + PDF + SVG):
  - king_r0_plot.{png,pdf,svg}
  - downsampling_plot.{png,pdf,svg}   (if downsampling enabled)
  - qc_panels.{png,pdf,svg}           (new; combines 3 QC panels)

QC panels include:
  A) Informative sites vs ANGSD parameter set (standard vs relaxed attempts)
  B) Per-sample depth distribution at informative sites (samtools depth; lightweight)
  C) Allele-balance (mixture heuristic) histograms from counts-based MAF (if available)

Notes:
  - PDF/SVG output is intended to be Illustrator-friendly.
"""

import argparse
import gzip
import os
import subprocess
import pathlib
import re
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Illustrator-friendly defaults
plt.rcParams["pdf.fonttype"] = 42      # embed TrueType
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"  # keep text as text


# ----------------------------
# Loaders
# ----------------------------
def load_ngsrelate_results(outdir: str) -> Optional[pd.DataFrame]:
    fn = os.path.join(outdir, "ngsrelate", "relatedness_results.tsv")
    return pd.read_csv(fn, sep="\t") if os.path.exists(fn) else None


def load_sample_stats(outdir: str) -> Optional[pd.DataFrame]:
    fn = os.path.join(outdir, "validation", "sample_stats.tsv")
    return pd.read_csv(fn, sep="\t") if os.path.exists(fn) else None


def load_downsampling_results(outdir: str) -> Optional[pd.DataFrame]:
    fn = os.path.join(outdir, "validation", "downsampling", "downsampling_results.tsv")
    return pd.read_csv(fn, sep="\t") if os.path.exists(fn) else None


def load_mixture_summary(outdir: str) -> Optional[str]:
    fn = os.path.join(outdir, "validation", "mixture", "mixture_summary.txt")
    if not os.path.exists(fn):
        return None
    with open(fn, "r") as f:
        return f.read()


def read_text_file(path: str, default: str = "") -> str:
    if not os.path.exists(path):
        return default
    with open(path, "r") as f:
        return f.read().strip()


def _pair_label(row: pd.Series) -> str:
    """Best-effort pair label across different column conventions."""
    a = None
    b = None
    for key in ["individual_A", "sample_1", "sample1", "ind1", "id1", "a", "bam1", "name1"]:
        if key in row and pd.notna(row.get(key)):
            a = str(row.get(key))
            break
    for key in ["individual_B", "sample_2", "sample2", "ind2", "id2", "b", "bam2", "name2"]:
        if key in row and pd.notna(row.get(key)):
            b = str(row.get(key))
            break
    if a is None:
        a = "?"
    if b is None:
        b = "?"
    return f"{a} vs {b}"


# ----------------------------
# Helpers
# ----------------------------
def ensure_dirs(outdir: str) -> None:
    os.makedirs(os.path.join(outdir, "reports"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "plots"), exist_ok=True)


def parse_plot_formats() -> List[str]:
    # Default: all three
    fmts = os.environ.get("PLOT_FORMATS", "png,pdf,svg")
    out = [x.strip().lower() for x in fmts.split(",") if x.strip()]
    # enforce uniqueness and safe ordering
    seen = set()
    ordered = []
    for f in ["png", "pdf", "svg"]:
        if f in out and f not in seen:
            ordered.append(f); seen.add(f)
    # include any extras after
    for f in out:
        if f not in seen:
            ordered.append(f); seen.add(f)
    return ordered


def _sanitize_svg_fonts(svg_path: str) -> None:
    """Replace common fallback fonts with Arial to avoid Illustrator import warnings."""
    repl = {
        "DejaVu Sans": "Arial",
        "Bitstream Vera Sans": "Arial",
        "Computer Modern Sans Serif": "Arial",
        "Geneva": "Arial",
        "Lucid": "Arial",
        "Avant Garde": "Arial",
        "Helvetica": "Arial",
        "Liberation Sans": "Arial",
    }
    p = pathlib.Path(svg_path)
    try:
        s = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return
    for k, v in repl.items():
        s = s.replace(k, v)
    # Normalize any remaining stacks that end in sans-serif
    s = re.sub(r"font-family:\s*'[^']*',\s*sans-serif", "font-family: Arial, sans-serif", s)
    s = re.sub(r"font-family:\s*[^;]*sans-serif", "font-family: Arial, sans-serif", s)
    p.write_text(s, encoding="utf-8")

def save_figure(fig: plt.Figure, basepath: str, formats: List[str], dpi_png: int = 200) -> List[str]:
    """Save figure to multiple formats (png/pdf/svg) with Illustrator-friendly defaults."""
    paths = []
    for ext in formats:
        out = f"{basepath}.{ext}"

        # Per-format typography (keep SVG text editable; embed TrueType in PDF)
        if ext in ("pdf", "svg"):
            plt.rcParams["pdf.fonttype"] = 42
            plt.rcParams["ps.fonttype"] = 42
            plt.rcParams["svg.fonttype"] = "none"

        if ext == "png":
            fig.savefig(out, dpi=dpi_png, bbox_inches="tight")
        else:
            fig.savefig(out, bbox_inches="tight")

        if ext == "svg":
            _sanitize_svg_fonts(out)

        paths.append(out)
    return paths




def stream_mafs_positions(mafs_gz: str, max_sites: int = 200_000) -> List[Tuple[str, int]]:
    """Read (chromo, pos) from mafs.gz without loading full file."""
    out: List[Tuple[str, int]] = []
    with gzip.open(mafs_gz, "rt") as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            out.append((chrom, pos))
            if len(out) >= max_sites:
                break
    return out


def write_bed_from_positions(positions: List[Tuple[str, int]], bed_path: str) -> None:
    with open(bed_path, "w") as out:
        for chrom, pos in positions:
            out.write(f"{chrom}\t{pos-1}\t{pos}\n")


def run_samtools_depth(bam: str, bed: str, mapq: int = 30, bq: int = 20) -> Optional[np.ndarray]:
    """
    Returns depth values for all positions in BED using samtools depth -a -b.
    If samtools fails/not installed, returns None.
    """
    cmd = ["samtools", "depth", "-a", "-d", "0", "-q", str(mapq), "-Q", str(bq), "-b", bed, bam]
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception:
        return None

    depths = []
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                depths.append(int(parts[2]))
            except ValueError:
                continue
    return np.array(depths, dtype=int)


def load_angsd_site_counts(outdir: str) -> Dict[str, int]:
    """
    Returns site counts for standard/relaxed attempts if present.
    Looks for:
      angsd/joint_gl.nsites
      angsd/joint_gl_relaxed1.nsites
      angsd/joint_gl_relaxed2.nsites
    """
    angsd_dir = os.path.join(outdir, "angsd")
    candidates = [
        ("standard", os.path.join(angsd_dir, "joint_gl.nsites")),
        ("relaxed1", os.path.join(angsd_dir, "joint_gl_relaxed1.nsites")),
        ("relaxed2", os.path.join(angsd_dir, "joint_gl_relaxed2.nsites")),
    ]
    out = {}
    for label, path in candidates:
        if os.path.exists(path):
            try:
                with open(path) as f:
                    out[label] = int(f.read().strip())
            except Exception:
                pass
    return out


def load_mixture_maf_values(outdir: str) -> Dict[str, pd.DataFrame]:
    """
    Loads counts-based MAF values per sample if mixture_analysis.py ran.
    Expects validation/mixture/<sample>/maf_values.tsv.gz with a 'maf' column.
    If present, also reads 'depth' or 'depth_major_minor' for depth-aware plots.
    """
    mix_dir = os.path.join(outdir, "validation", "mixture")
    out: Dict[str, pd.DataFrame] = {}
    if not os.path.isdir(mix_dir):
        return out
    for sample in os.listdir(mix_dir):
        p = os.path.join(mix_dir, sample, "maf_values.tsv.gz")
        if not os.path.exists(p):
            continue
        try:
            df = pd.read_csv(p, sep="\t")
            if "maf" not in df.columns or len(df) == 0:
                continue
            if "depth" not in df.columns and "depth_major_minor" in df.columns:
                df = df.rename(columns={"depth_major_minor": "depth"})
            keep = ["maf"]
            if "depth" in df.columns:
                keep.append("depth")
            out[sample] = df[keep].dropna(subset=["maf"]).copy()
        except Exception:
            continue
    return out






# ----------------------------
# Plots (publication-ready)
# ----------------------------
def plot_king_r0(results_df: pd.DataFrame, outdir: str, formats: List[str]) -> None:
    base = os.path.join(outdir, "plots", "king_r0_plot")

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    categories = [
        (0.354, 0.6, 'Identical', '#2ecc71', 0.15),
        (0.177, 0.354, '1st degree', '#3498db', 0.15),
        (0.0884, 0.177, '2nd degree', '#f39c12', 0.15),
        (0.0442, 0.0884, '3rd degree', '#e74c3c', 0.10),
        (-0.1, 0.0442, 'Unrelated', '#95a5a6', 0.08),
    ]
    for ymin, ymax, label, color, alpha in categories:
        ax.axhspan(ymin, ymax, alpha=alpha, color=color, label=label)
        ax.axhline(y=ymax, color=color, linestyle='--', alpha=0.4, linewidth=0.8)

    for _, row in results_df.iterrows():
        king = row.get('KING_kinship', np.nan)
        r0 = row.get('R0', np.nan)
        label = _pair_label(row)
        if pd.notna(r0) and pd.notna(king):
            ax.scatter(r0, king, s=120, c='black', zorder=5, edgecolors='white', linewidth=1.5)
            ax.annotate(label, (r0, king), textcoords="offset points",
                        xytext=(10, 5), fontsize=8, ha='left')
        elif pd.notna(king):
            ax.scatter(0, king, s=120, c='black', zorder=5, marker='D',
                       edgecolors='white', linewidth=1.5)
            ax.annotate(label + "\n(R0 unavailable)", (0, king),
                        textcoords="offset points", xytext=(10, 5), fontsize=8, ha='left')

    ax.set_xlabel('R0 (fraction of homozygous-opposite sites)', fontsize=12)
    ax.set_ylabel('KING-robust kinship coefficient', fontsize=12)
    ax.set_title('Pairwise Relatedness Estimates\n(ANGSD + NgsRelate, transversions only)',
                 fontsize=13)
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

    ax.set_xlim(-0.05, max(1.05, ax.get_xlim()[1]))
    king_vals = results_df['KING_kinship'].dropna() if 'KING_kinship' in results_df.columns else pd.Series(dtype=float)
    if not king_vals.empty:
        ymin_plot = min(-0.1, float(king_vals.min()) - 0.05)
        ymax_plot = max(0.6, float(king_vals.max()) + 0.05)
        ax.set_ylim(ymin_plot, ymax_plot)
    else:
        ax.set_ylim(-0.1, 0.6)

    plt.tight_layout()
    save_figure(fig, base, formats=formats, dpi_png=200)
    plt.close(fig)


def plot_downsampling(ds_df: pd.DataFrame, outdir: str, formats: List[str]) -> None:
    base = os.path.join(outdir, "plots", "downsampling_plot")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ds_valid = ds_df[ds_df['king'] != 'NA'].copy()
    if ds_valid.empty:
        ax1.text(0.5, 0.5, 'No valid downsampling data', transform=ax1.transAxes,
                 ha='center', va='center')
        ax2.text(0.5, 0.5, 'No valid downsampling data', transform=ax2.transAxes,
                 ha='center', va='center')
    else:
        ds_valid['king'] = ds_valid['king'].astype(float)
        ds_valid['fraction'] = ds_valid['fraction'].astype(float)
        ds_valid['n_sites'] = ds_valid['n_sites'].astype(int)

        ax1.plot(ds_valid['fraction'], ds_valid['king'], 'ko-', markersize=8)
        for _, row in ds_valid.iterrows():
            ax1.annotate(row['label'], (row['fraction'], row['king']),
                         textcoords="offset points", xytext=(5, 5), fontsize=8)

        for thresh, label, color in [(0.354, 'Identical', '#2ecc71'),
                                     (0.177, '1st degree', '#3498db'),
                                     (0.0884, '2nd degree', '#f39c12'),
                                     (0.0442, '3rd degree', '#e74c3c')]:
            ax1.axhline(y=thresh, color=color, linestyle='--', alpha=0.5, linewidth=0.8,
                        label=f'{label} ({thresh})')

        ax1.set_xlabel('Downsampling Fraction', fontsize=11)
        ax1.set_ylabel('KING-robust kinship', fontsize=11)
        ax1.set_title('KING Kinship vs Coverage', fontsize=12)
        ax1.legend(fontsize=8, loc='best')

        ax2.plot(ds_valid['fraction'], ds_valid['n_sites'], 'bs-', markersize=8)
        for _, row in ds_valid.iterrows():
            ax2.annotate(f"{row['n_sites']}", (row['fraction'], row['n_sites']),
                         textcoords="offset points", xytext=(5, 5), fontsize=8)

        ax2.set_xlabel('Downsampling Fraction', fontsize=11)
        ax2.set_ylabel('Number of Informative Sites', fontsize=11)
        ax2.set_title('Informative Sites vs Coverage', fontsize=12)

    plt.tight_layout()
    save_figure(fig, base, formats=formats, dpi_png=200)
    plt.close(fig)


def plot_qc_panels(outdir: str, bams_file: str, formats: List[str]) -> None:
    """
    Generates a 4-panel QC figure (2x2):
      A) Sites found by ANGSD parameter set
      B) Depth at informative sites per sample (samtools depth on informative sites)
      C) Allele-balance histograms (minor/(major+minor)) from mixture step
      D) Allele-balance vs depth (helps distinguish low-depth noise vs mixture/reference bias)
    """
    base = os.path.join(outdir, "plots", "qc_panels")
    angsd_dir = os.path.join(outdir, "angsd")
    mafs_gz = os.path.join(angsd_dir, "joint_gl.mafs.gz")
    chosen = read_text_file(os.path.join(angsd_dir, "param_set_used.txt"), default="standard")
    site_counts = load_angsd_site_counts(outdir)

    stats_df = load_sample_stats(outdir)
    if stats_df is None:
        return

    # Prepare informative sites BED
    bed_path = os.path.join(outdir, "plots", "informative_sites.bed")
    depths_by_sample: Dict[str, Optional[np.ndarray]] = {}

    positions = []
    if os.path.exists(mafs_gz):
        positions = stream_mafs_positions(mafs_gz, max_sites=200_000)
        if positions:
            write_bed_from_positions(positions, bed_path)

    # Load BAM list
    bam_paths = []
    with open(bams_file, "r") as f:
        for line in f:
            p = line.strip()
            if p:
                bam_paths.append(p)

    # Compute depth distributions (only at informative sites)
    if positions:
        for bam in bam_paths:
            sample = os.path.basename(bam).replace(".bam", "")
            depths = run_samtools_depth(bam, bed_path, mapq=30, bq=20)
            depths_by_sample[sample] = depths
    else:
        for bam in bam_paths:
            sample = os.path.basename(bam).replace(".bam", "")
            depths_by_sample[sample] = None

    # Mixture MAF tables (counts-based)
    maf_tables = load_mixture_maf_values(outdir)

    fig, axes = plt.subplots(2, 2, figsize=(16, 8))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # Panel A
    labels = []
    values = []
    for key in ["standard", "relaxed1", "relaxed2"]:
        if key in site_counts:
            labels.append(key)
            values.append(site_counts[key])

    if not labels:
        axA.text(0.5, 0.5, "No ANGSD site count info", ha="center", va="center", transform=axA.transAxes)
    else:
        x = np.arange(len(labels))
        axA.bar(x, values)
        axA.set_xticks(x)
        axA.set_xticklabels(labels, rotation=0)
        for i, v in enumerate(values):
            axA.text(i, v + max(1, 0.02 * max(values)), str(v), ha="center", va="bottom", fontsize=9)
        axA.set_ylabel("Informative sites (n)")
        axA.set_title("A) Sites vs ANGSD parameter set")
        axA.text(0.02, 0.92, f"chosen: {chosen}", transform=axA.transAxes, fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.8", alpha=0.9))

    # Panel B
    axB.set_title("B) Depth at informative sites")
    axB.set_xlabel("Depth (samtools depth)")
    axB.set_ylabel("Sites")

    have_depth = False
    for sample, depths in depths_by_sample.items():
        if depths is None or len(depths) == 0:
            continue
        have_depth = True
        max_d = int(np.max(depths))
        bins = np.arange(0, max(6, max_d + 2)) - 0.5
        axB.hist(depths, bins=bins, histtype="step", linewidth=2, label=sample)

    if not have_depth:
        axB.text(0.5, 0.5, "Depth unavailable\n(samtools missing or no sites)", ha="center",
                 va="center", transform=axB.transAxes)
    else:
        axB.legend(fontsize=8, framealpha=0.9)

    # Panel C
    axC.set_title("C) Allele-balance (minor/(major+minor))")
    axC.set_xlabel("Minor allele fraction")
    axC.set_ylabel("Sites")

    have_maf = False
    for sample, df in maf_tables.items():
        if df is None or len(df) == 0 or "maf" not in df.columns:
            continue
        have_maf = True
        axC.hist(df["maf"].to_numpy(dtype=float), bins=25, range=(0, 0.5),
                 histtype="step", linewidth=2, label=sample)

    if not have_maf:
        axC.text(0.5, 0.5, "No mixture outputs\n(or insufficient data)", ha="center",
                 va="center", transform=axC.transAxes)
    else:
        axC.set_xlim(0, 0.5)
        axC.legend(fontsize=8, framealpha=0.9)

    # Panel D
    axD.set_title("D) Allele-balance vs depth")
    axD.set_xlabel("Depth (major+minor)")
    axD.set_ylabel("Minor allele fraction")

    have_scatter = False
    for sample, df in maf_tables.items():
        if df is None or len(df) == 0:
            continue
        if "maf" not in df.columns or "depth" not in df.columns:
            continue
        have_scatter = True
        axD.scatter(df["depth"].to_numpy(dtype=float), df["maf"].to_numpy(dtype=float),
                    s=12, alpha=0.35, label=sample, edgecolors="none")

    if not have_scatter:
        axD.text(0.5, 0.5, "No depth+MAF table\n(mixture step low coverage?)",
                 ha="center", va="center", transform=axD.transAxes)
    else:
        axD.set_ylim(-0.02, 0.52)
        axD.set_xlim(left=0)
        axD.legend(fontsize=8, framealpha=0.9, loc="lower right")

    plt.suptitle("Relatedness QC panels", fontsize=14, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_figure(fig, base, formats=formats, dpi_png=200)
    plt.close(fig)




# ----------------------------
# Interpretation
# ----------------------------
def generate_recommendation(results_df: Optional[pd.DataFrame]) -> Tuple[str, str]:
    if results_df is None or results_df.empty:
        return "INCONCLUSIVE", "No relatedness results available."

    row = results_df.iloc[0]
    category = row.get('category', 'UNKNOWN')
    confidence = str(row.get('confidence', 'UNKNOWN'))
    king = row.get('KING_kinship', float('nan'))
    n_sites = int(row.get('n_sites', 0)) if str(row.get('n_sites', '')).strip() != "" else 0

    if 'VERY_LOW' in confidence:
        rec = "INCONCLUSIVE"
        expl = (
            f"With only {n_sites} informative sites, the data is insufficient for a reliable determination. "
            "The KING estimate may be highly variable, especially under cross-species mapping."
        )
    elif category == "IDENTICAL":
        rec = "SAME INDIVIDUAL"
        expl = f"KING kinship = {king:.4f} (> 0.354 identity threshold)."
    elif category == "UNRELATED":
        rec = "DIFFERENT INDIVIDUALS"
        expl = f"KING kinship = {king:.4f} (< 0.0442 unrelated threshold)."
    elif category in ("FIRST_DEGREE", "SECOND_DEGREE", "THIRD_DEGREE"):
        rec = "DIFFERENT INDIVIDUALS (possibly related)"
        expl = f"KING kinship = {king:.4f} indicates {category.lower().replace('_', ' ')} relatedness."
    else:
        rec = "INCONCLUSIVE"
        expl = f"Unable to classify relationship (KING = {king})."

    if 'LOW' in confidence and 'VERY' not in confidence:
        expl += f" NOTE: Low confidence due to limited informative sites ({n_sites})."

    return rec, expl


def main():
    ap = argparse.ArgumentParser(description="Generate summary report and plots")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bams", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--n-samples", type=int, required=True)
    ap.add_argument("--downsample", default="false")
    ap.add_argument("--mixture-check", default="false")
    args = ap.parse_args()

    ensure_dirs(args.outdir)
    formats = parse_plot_formats()

    results_df = load_ngsrelate_results(args.outdir)
    stats_df = load_sample_stats(args.outdir)
    ds_df = load_downsampling_results(args.outdir) if args.downsample.lower() == "true" else None
    mixture_text = load_mixture_summary(args.outdir) if args.mixture_check.lower() == "true" else None

    # Metadata
    angsd_dir = os.path.join(args.outdir, "angsd")
    param_set = read_text_file(os.path.join(angsd_dir, "param_set_used.txt"), default="standard")
    nsites_file = os.path.join(angsd_dir, "joint_gl.nsites")
    try:
        n_sites = int(read_text_file(nsites_file, default="0") or "0")
    except Exception:
        n_sites = 0

    # Plots
    if results_df is not None and not results_df.empty:
        plot_king_r0(results_df, args.outdir, formats=formats)

    if ds_df is not None and not ds_df.empty:
        plot_downsampling(ds_df, args.outdir, formats=formats)

    # New QC panels (always attempt; handles missing components gracefully)
    plot_qc_panels(args.outdir, args.bams, formats=formats)

    # Recommendation
    recommendation, explanation = generate_recommendation(results_df)

    # Report
    report_path = os.path.join(args.outdir, "reports", "summary_report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 72 + "\n")
        f.write("  PalaeoKin aDNA RELATEDNESS PIPELINE — SUMMARY REPORT\n")
        f.write("=" * 72 + "\n")
        f.write(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 72 + "\n\n")

        f.write("┌" + "─" * 70 + "┐\n")
        f.write(f"│  RECOMMENDATION: {recommendation:<52}│\n")
        f.write("└" + "─" * 70 + "┘\n\n")
        f.write(f"  {explanation}\n\n")

        f.write("─" * 72 + "\n")
        f.write("  SAMPLE INFORMATION\n")
        f.write("─" * 72 + "\n")
        if stats_df is not None:
            for _, row in stats_df.iterrows():
                f.write(f"  {row['sample']}\n")
                f.write(f"    Reads:    {int(row['total_reads']):>12,}\n")
                f.write(f"    Avg len:  {row['avg_read_length']:>12} bp\n")
                f.write(f"    Coverage: {row['est_coverage']:>12.6f}x\n")
                f.write(f"    BAM:      {row['bam_path']}\n\n")
        f.write(f"  Reference: {args.ref}\n\n")

        f.write("─" * 72 + "\n")
        f.write("  ANGSD ANALYSIS\n")
        f.write("─" * 72 + "\n")
        f.write(f"  Parameter set used:  {param_set}\n")
        f.write(f"  Informative sites:   {n_sites}\n\n")

        f.write("─" * 72 + "\n")
        f.write("  RELATEDNESS RESULTS (NgsRelate)\n")
        f.write("─" * 72 + "\n")
        if results_df is None or results_df.empty:
            f.write("  No relatedness results available.\n\n")
        else:
            row = results_df.iloc[0]
            f.write(f"  {_pair_label(row)}\n")
            f.write(f"    KING kinship:    {row.get('KING_kinship','NA')}\n")
            f.write(f"    R0:              {row.get('R0','NA')}\n")
            f.write(f"    R1:              {row.get('R1','NA')}\n")
            f.write(f"    Theta:           {row.get('theta','NA')}\n")
            f.write(f"    Sites used:      {row.get('n_sites','NA')}\n")
            f.write(f"    Classification:  {row.get('category','NA')}\n")
            f.write(f"    Confidence:      {row.get('confidence','NA')}\n\n")

        f.write("─" * 72 + "\n")
        f.write("  INTERPRETATION GUIDE\n")
        f.write("─" * 72 + "\n")
        f.write("  KING-robust kinship (rule-of-thumb thresholds):\n")
        f.write("    ~0.50  : identical / duplicate library / same individual (needs many sites)\n")
        f.write("    ~0.25  : 1st degree relatives (parent-offspring or full siblings)\n")
        f.write("    ~0.125 : 2nd degree relatives\n")
        f.write("    ~0.0625: 3rd degree relatives\n")
        f.write("    ~0.00  : unrelated\n\n")
        f.write("  R0 / R1:\n")
        f.write("    R0 = fraction of homozygous-opposite sites (close relatives tend to have low R0).\n")
        f.write("    R1 = allele-sharing statistic used with R0 + KING.\n\n")
        f.write("  Sites used (n_sites) drives confidence:\n")
        f.write("    Very low n_sites (e.g., <~100) makes estimates unstable; interpret conservatively.\n")
        f.write("    Cross-species mapping can bias all statistics.\n\n")


        if ds_df is not None and not ds_df.empty:
            f.write("─" * 72 + "\n")
            f.write("  DOWNSAMPLING VALIDATION\n")
            f.write("─" * 72 + "\n")
            f.write("  See plots/downsampling_plot.* for stability across downsampling fractions.\n\n")

        if mixture_text:
            f.write("─" * 72 + "\n")
            f.write("  MIXTURE / ALLELE-BALANCE CHECK\n")
            f.write("─" * 72 + "\n")
            f.write(mixture_text.strip() + "\n\n")

        f.write("─" * 72 + "\n")
        f.write("  PLOTS\n")
        f.write("─" * 72 + "\n")
        f.write("  - plots/king_r0_plot.{png,pdf,svg}\n")
        f.write("  - plots/downsampling_plot.{png,pdf,svg} (if enabled)\n")
        f.write("  - plots/qc_panels.{png,pdf,svg}\n")

    print(f"Summary report: {report_path}")
    print(f"Plots: {os.path.join(args.outdir, 'plots')}")


if __name__ == "__main__":
    main()
