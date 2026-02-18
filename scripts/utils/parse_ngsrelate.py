#!/usr/bin/env python3
"""
parse_ngsrelate.py — Parse NgsRelate output and classify relationships.

Reads NgsRelate output TSV and classifies each pairwise comparison using:
  - KING-robust kinship coefficient (rab)
  - R0, R1 statistics to distinguish relationship types
  - Confidence rating based on number of informative sites

KING-robust thresholds (Manichaikul et al. 2010; Waples et al. 2019):
  > 0.354      : Identical / monozygotic twins
  0.177 - 0.354: First-degree relatives (parent-offspring or full siblings)
  0.0884- 0.177: Second-degree relatives (half-siblings, avuncular, grandparent)
  0.0442-0.0884: Third-degree relatives
  < 0.0442     : Unrelated

R0/R1 distinguish within first-degree:
  - Parent-offspring: R0 ≈ 0, R1 high
  - Full siblings: R0 > 0, R1 moderate
"""

import argparse
import os
import sys

import pandas as pd


def parse_ngsrelate_output(input_file):
    """Parse NgsRelate output TSV."""
    df = pd.read_csv(input_file, sep='\t')

    # NgsRelate column names vary; identify key columns
    # Standard columns: a, b, rab, Fa, Fb, theta, inbreed_a, inbreed_b, R0, R1, KING, nSNPs
    # Some versions use different headers

    # Map common column name variants
    col_map = {}
    for col in df.columns:
        col_lower = col.lower().strip()
        if col_lower in ('rab', 'king', 'king_robust'):
            col_map['king'] = col
        elif col_lower == 'r0':
            col_map['R0'] = col
        elif col_lower == 'r1':
            col_map['R1'] = col
        elif col_lower in ('nsnps', 'nsnp', 'nsites', 'nsit'):
            col_map['nSNPs'] = col
        elif col_lower == 'a':
            col_map['a'] = col
        elif col_lower == 'b':
            col_map['b'] = col
        elif col_lower == 'ida':
            col_map['a'] = col
        elif col_lower == 'idb':
            col_map['b'] = col
        elif col_lower == 'theta':
            col_map['theta'] = col

    return df, col_map


def classify_relationship(king, r0=None, r1=None):
    """Classify relationship based on KING-robust kinship."""
    if king > 0.354:
        category = "IDENTICAL"
        detail = "Same individual or monozygotic twins"
    elif king > 0.177:
        category = "FIRST_DEGREE"
        if r0 is not None and r1 is not None:
            if r0 < 0.1:
                detail = "Parent-offspring (R0 ≈ 0)"
            else:
                detail = "Full siblings (R0 > 0)"
        else:
            detail = "First-degree relatives"
    elif king > 0.0884:
        category = "SECOND_DEGREE"
        detail = "Second-degree relatives (half-sibs, avuncular, grandparent)"
    elif king > 0.0442:
        category = "THIRD_DEGREE"
        detail = "Third-degree relatives"
    else:
        category = "UNRELATED"
        detail = "Unrelated individuals"
    return category, detail


def confidence_rating(n_sites, param_set="standard"):
    """Rate confidence based on number of informative sites."""
    if n_sites < 100:
        rating = "VERY_LOW"
        note = "< 100 informative sites; results are unreliable"
    elif n_sites < 500:
        rating = "LOW"
        note = "100-500 informative sites; results should be interpreted with caution"
    elif n_sites < 2000:
        rating = "MODERATE"
        note = "500-2000 informative sites; results are reasonably reliable"
    else:
        rating = "HIGH"
        note = ">= 2000 informative sites; results are reliable"

    if param_set != "standard":
        rating = rating + " (RELAXED_PARAMS)"
        note += f"; used relaxed parameters ({param_set})"

    return rating, note


def main():
    parser = argparse.ArgumentParser(description="Parse NgsRelate output")
    parser.add_argument("--input", required=True, help="NgsRelate output TSV")
    parser.add_argument("--n-sites", type=int, required=True,
                        help="Number of ANGSD sites used")
    parser.add_argument("--param-set", default="standard",
                        help="Parameter set used (standard/relaxed_pval_maf/relaxed_all)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    df, col_map = parse_ngsrelate_output(args.input)

    if df.empty:
        print("ERROR: NgsRelate output is empty", file=sys.stderr)
        sys.exit(1)

    # Build results
    results = []
    for _, row in df.iterrows():
        # Get individual IDs
        ind_a = row.get(col_map.get('a', 'a'), 'Unknown_A')
        ind_b = row.get(col_map.get('b', 'b'), 'Unknown_B')

        # Get KING kinship
        king_col = col_map.get('king')
        if king_col and king_col in row:
            king = float(row[king_col])
        else:
            # Fall back to column index 6 (0-based) which is typically rab/KING
            king = float(row.iloc[6]) if len(row) > 6 else float('nan')

        # Get R0, R1
        r0 = None
        r1 = None
        r0_col = col_map.get('R0')
        r1_col = col_map.get('R1')
        if r0_col and r0_col in row:
            try:
                r0 = float(row[r0_col])
            except (ValueError, TypeError):
                pass
        if r1_col and r1_col in row:
            try:
                r1 = float(row[r1_col])
            except (ValueError, TypeError):
                pass

        # Get nSNPs from NgsRelate (may differ from ANGSD site count)
        nsnps_col = col_map.get('nSNPs')
        if nsnps_col and nsnps_col in row:
            try:
                nsnps = int(row[nsnps_col])
            except (ValueError, TypeError):
                nsnps = args.n_sites
        else:
            nsnps = args.n_sites

        # Classify
        category, detail = classify_relationship(king, r0, r1)
        conf_rating, conf_note = confidence_rating(nsnps, args.param_set)

        # Get theta if available
        theta_col = col_map.get('theta')
        theta = None
        if theta_col and theta_col in row:
            try:
                theta = float(row[theta_col])
            except (ValueError, TypeError):
                pass

        results.append({
            'individual_A': ind_a,
            'individual_B': ind_b,
            'KING_kinship': king,
            'R0': r0,
            'R1': r1,
            'theta': theta,
            'n_sites': nsnps,
            'category': category,
            'detail': detail,
            'confidence': conf_rating,
            'confidence_note': conf_note,
        })

    # Write results TSV
    results_df = pd.DataFrame(results)
    results_tsv = os.path.join(args.outdir, "relatedness_results.tsv")
    results_df.to_csv(results_tsv, sep='\t', index=False)
    print(f"Results written to: {results_tsv}")

    # Write human-readable report
    report_file = os.path.join(args.outdir, "relatedness_report.txt")
    with open(report_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("NgsRelate Relatedness Report\n")
        f.write("=" * 70 + "\n\n")

        for r in results:
            f.write(f"Pair: {r['individual_A']}  vs  {r['individual_B']}\n")
            f.write("-" * 50 + "\n")
            f.write(f"  KING-robust kinship:  {r['KING_kinship']:.6f}\n")
            if r['R0'] is not None:
                f.write(f"  R0:                   {r['R0']:.6f}\n")
            if r['R1'] is not None:
                f.write(f"  R1:                   {r['R1']:.6f}\n")
            if r['theta'] is not None:
                f.write(f"  Theta:                {r['theta']:.6f}\n")
            f.write(f"  Informative sites:    {r['n_sites']}\n")
            f.write(f"\n")
            f.write(f"  Classification:       {r['category']}\n")
            f.write(f"  Detail:               {r['detail']}\n")
            f.write(f"  Confidence:           {r['confidence']}\n")
            f.write(f"  Note:                 {r['confidence_note']}\n")
            f.write("\n")

        f.write("=" * 70 + "\n")
        f.write("KING-robust kinship thresholds:\n")
        f.write("  > 0.354       : Identical / same individual\n")
        f.write("  0.177 - 0.354 : First-degree (parent-offspring / full siblings)\n")
        f.write("  0.0884 - 0.177: Second-degree (half-siblings, avuncular)\n")
        f.write("  0.0442 - 0.0884: Third-degree relatives\n")
        f.write("  < 0.0442      : Unrelated\n")
        f.write("=" * 70 + "\n")

    print(f"Report written to: {report_file}")

    # Print summary to stdout
    for r in results:
        print(f"\n  {r['individual_A']} vs {r['individual_B']}: "
              f"KING={r['KING_kinship']:.4f} -> {r['category']} "
              f"({r['confidence']}, {r['n_sites']} sites)")


if __name__ == "__main__":
    main()
