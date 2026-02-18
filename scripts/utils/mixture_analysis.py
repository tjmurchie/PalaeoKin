#!/usr/bin/env python3
"""Compute simple allele-balance / mixture heuristics from ANGSD dumpCounts.

Input
-----
--counts : path to ANGSD -dumpCounts output (optionally .gz)

Output (written to --outdir)
---------------------------
maf_values.tsv.gz :
  chr  pos  depth  major  minor  maf
    - depth = total depth (A+C+G+T)
    - major/minor = bases (A/C/G/T) with highest/second-highest observed counts
    - maf = minor/(major+minor)  (0..0.5), ignoring 3rd/4th alleles by renormalising

mixture_report.txt :
  summary stats useful for quick interpretation

maf_histogram.* :
  histogram plot (formats configurable)

Notes
-----
ANGSD dumpCounts formats vary slightly depending on flags/version. This script
tries to robustly locate the final A/C/G/T count columns by defaulting to the
*last four* numeric columns on each line (works for common layouts like:
  chr pos A C G T
  chr pos depth A C G T
  chr pos ... A C G T
).
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

# Matplotlib is only used for the quick histogram output (optional but handy).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASES = np.array(["A", "C", "G", "T"], dtype=object)


def _open_text_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def _parse_ints(tokens: List[str]) -> Optional[List[int]]:
    try:
        return [int(x) for x in tokens]
    except Exception:
        return None


def _extract_acgt_counts(parts: List[str]) -> Optional[List[int]]:
    """Return [A,C,G,T] counts from a split line, or None if not parseable.

    Supports:
      - ANGSD dumpCounts style lines (chr pos ... A C G T)
      - ANGSD .counts.gz style lines with ONLY A C G T columns (common with -doCounts 1)
    """
    if len(parts) < 4:
        return None

    # Most robust: last 4 columns are A/C/G/T counts.
    cand = _parse_ints(parts[-4:])
    if cand is not None and len(cand) == 4:
        return cand

    # Fallbacks for older / unusual layouts.
    cand = _parse_ints(parts[2:6])
    if cand is not None and len(cand) == 4:
        return cand

    if len(parts) >= 7:
        cand = _parse_ints(parts[3:7])  # e.g., chr pos depth A C G T
        if cand is not None and len(cand) == 4:
            return cand

    return None


def _savefig_multi(fig, outprefix: str, formats: str):
    fmt_list = [f.strip().lower() for f in formats.split(",") if f.strip()]
    if not fmt_list:
        fmt_list = ["png"]
    for fmt in fmt_list:
        fig.savefig(f"{outprefix}.{fmt}", dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="ANGSD dumpCounts output OR A/C/G/T-only .counts(.gz)")
    ap.add_argument("--pos", default=None, help="Optional positions file (chr	pos). If not provided and counts lacks chr/pos, will try sibling .pos(.gz).")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--sample", default=None, help="Sample name (optional; used in report)")
    ap.add_argument(
        "--formats",
        default=None,
        help="Plot formats, e.g. 'png,pdf,svg'. If not set, uses $PLOT_FORMATS or defaults to png.",
    )
    ap.add_argument(
        "--min-depth",
        type=int,
        default=2,
        help="Minimum total depth (A+C+G+T) to retain a site in maf_values.tsv.gz (default: 2)",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Allow --formats to override the env var (the pipeline previously passed --formats).
    plot_formats = args.formats or os.environ.get("PLOT_FORMATS", "png")

    rows: List[Tuple[str, int, int, str, str, float, int, int]] = []
    n_read = 0
    n_skipped_parse = 0
    n_skipped_depth = 0

    # Helper: if counts file is A/C/G/T-only, pull chr/pos from a separate .pos file.
    pos_path = args.pos
    if pos_path is None:
        # Try sibling file: <prefix>.counts(.gz) -> <prefix>.pos(.gz)
        p = str(args.counts)
        # Strip .gz for probing
        base = p[:-3] if p.endswith(".gz") else p
        if base.endswith(".counts"):
            cand = base[:-len(".counts")] + ".pos"
            if os.path.exists(cand):
                pos_path = cand
            elif os.path.exists(cand + ".gz"):
                pos_path = cand + ".gz"

    pos_iter = None
    pos_used = False
    if pos_path is not None and os.path.exists(pos_path):
        def _iter_pos(ppath: str):
            with _open_text_maybe_gz(ppath) as pf:
                for pline in pf:
                    pline = pline.strip()
                    if not pline or pline.startswith("#"):
                        continue
                    pparts = pline.split()
                    # skip header like: chr pos
                    if len(pparts) >= 2 and pparts[0].lower() == "chr" and pparts[1].lower() == "pos":
                        continue
                    if len(pparts) < 2:
                        continue
                    try:
                        yield (pparts[0], int(pparts[1]))
                    except Exception:
                        continue
        pos_iter = _iter_pos(pos_path)
        pos_used = True

    with _open_text_maybe_gz(args.counts) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()

            # Skip obvious headers (e.g., ind0_A ind0_C ind0_G ind0_T)
            if any(re.search(r"[A-Za-z_]", tok) for tok in parts) and _parse_ints(parts[-4:]) is None:
                continue

            acgt = _extract_acgt_counts(parts)
            if acgt is None:
                n_skipped_parse += 1
                continue

            # Determine whether chr/pos are embedded (dumpCounts) or must be pulled from .pos
            chrom = None
            pos = None
            if len(parts) >= 2:
                try:
                    # dumpCounts has chr + pos as first two columns
                    pos = int(parts[1])
                    chrom = parts[0]
                    # sanity: if the line is too short to be dumpCounts, ignore chrom/pos guess
                    if len(parts) < 6:
                        chrom = None
                        pos = None
                except Exception:
                    chrom = None
                    pos = None

            if chrom is None or pos is None:
                # counts-only mode
                if pos_iter is not None:
                    try:
                        chrom, pos = next(pos_iter)
                    except StopIteration:
                        chrom, pos = ("NA", n_read + 1)
                else:
                    chrom, pos = ("NA", n_read + 1)

            depth = int(sum(acgt))
            n_read += 1

            if depth < args.min_depth:
                n_skipped_depth += 1
                continue

            # Identify major/minor alleles by observed counts.
            order = np.argsort(acgt)[::-1]
            major_i, minor_i = int(order[0]), int(order[1])
            major_c, minor_c = int(acgt[major_i]), int(acgt[minor_i])

            denom = major_c + minor_c
            maf = float(minor_c / denom) if denom > 0 else float("nan")

            rows.append((chrom, pos, depth, str(BASES[major_i]), str(BASES[minor_i]), maf, major_c, minor_c))
    df = pd.DataFrame(
        rows,
        columns=["chr", "pos", "depth", "major", "minor", "maf", "major_count", "minor_count"],
    )
    out_table = outdir / "maf_values.tsv.gz"
    df.to_csv(out_table, sep="\t", index=False, compression="gzip")

    n_total = int(len(df))
    n_inform = int((df["minor_count"] > 0).sum())

    # Compute MAF summary primarily on informative sites (minor>0), but report totals too.
    df_inf = df[df["minor_count"] > 0].copy()
    n_inf = int(len(df_inf))
    mean_maf = float(df_inf["maf"].mean()) if n_inf else float("nan")
    med_maf = float(df_inf["maf"].median()) if n_inf else float("nan")
    std_maf = float(df_inf["maf"].std(ddof=0)) if n_inf else float("nan")
    frac_near = float((df_inf["maf"] >= 0.45).mean()) if n_inf else float("nan")
    frac_int = float(((df_inf["maf"] >= 0.05) & (df_inf["maf"] <= 0.40)).mean()) if n_inf else float("nan")

    rep: List[str] = []
    rep.append("Mixture / allele-balance diagnostic (counts-based)")
    rep.append(f"Counts file: {args.counts}")
    if pos_used: rep.append(f"Positions file: {pos_path}")
    else: rep.append("Positions: embedded in counts file or not available")
    if args.sample:
        rep.append(f"Sample: {args.sample}")
    rep.append(f"Counts file: {args.counts}")
    rep.append("")
    rep.append("Parsing / filtering:")
    rep.append(f"  Lines read: {n_read}")
    rep.append(f"  Lines skipped (unparseable): {n_skipped_parse}")
    rep.append(f"  Sites retained (depth>={args.min_depth}): {n_total}")
    rep.append(f"  Sites with minor>0 (informative): {n_inform}")
    rep.append("")
    rep.append("Summary (informative sites only):")
    rep.append(f"  Mean MAF: {mean_maf:.4f}")
    rep.append(f"  Median MAF: {med_maf:.4f}")
    rep.append(f"  Std MAF: {std_maf:.4f}")
    rep.append(f"  Fraction near 0.5 (>=0.45): {frac_near:.3%}")
    rep.append(f"  Fraction intermediate (0.05–0.40): {frac_int:.3%}")
    rep.append("")
    rep.append("How to interpret:")
    rep.append("  - maf = minor/(major+minor), 0–0.5. Single-individual informative sites often cluster near ~0.5 when depth is sufficient.")
    rep.append("  - At low depth, discrete values (1/3≈0.33, 2/5=0.40, etc.) are expected from sampling noise.")
    rep.append("  - Mixture tends to produce excess intermediate MAF at higher depth; see QC Panel D (MAF vs depth).")

    (outdir / "mixture_report.txt").write_text("\n".join(rep), encoding="utf-8")

    # Quick histogram (informative sites if any; otherwise all retained).
    plot_df = df_inf if n_inf else df
    if len(plot_df) > 0:
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax.hist(plot_df["maf"].astype(float), bins=np.linspace(0, 0.5, 51), histtype="step", linewidth=2)
        ax.set_title("Allele-balance histogram (counts-based)")
        ax.set_xlabel("Minor allele fraction")
        ax.set_ylabel("Sites")
        fig.tight_layout()
        _savefig_multi(fig, str(outdir / "maf_histogram"), plot_formats)
        plt.close(fig)

    print(
        f"[mixture_analysis.py] Wrote {out_table} "
        f"(sites_retained={n_total}, informative_minor>0={n_inform}, min_depth={args.min_depth})"
    )


if __name__ == "__main__":
    main()