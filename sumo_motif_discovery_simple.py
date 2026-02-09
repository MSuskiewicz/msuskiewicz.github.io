#!/usr/bin/env python3
"""
SUMO Site Analyzer - De Novo Motif Discovery (Simplified Version)
==================================================================
Simple version for discovering enriched sequence motifs around SUMO sites.

Usage:
    python sumo_motif_discovery_simple.py <core_output_excel> [output_file]
"""

import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
from collections import Counter

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

SCORE_VERY_HIGH_THRESHOLD = 600
SCORE_HIGH_THRESHOLD = 400
SCORE_MEDIUM_THRESHOLD = 200

MIN_MOTIF_COUNT = 5
FDR_THRESHOLD = 0.1

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

# Background amino acid frequencies (UniProt/SwissProt)
BACKGROUND_AA_FREQ = {
    'A': 0.0825, 'C': 0.0137, 'D': 0.0545, 'E': 0.0675,
    'F': 0.0386, 'G': 0.0707, 'H': 0.0227, 'I': 0.0596,
    'K': 0.0584, 'L': 0.0966, 'M': 0.0242, 'N': 0.0406,
    'P': 0.0470, 'Q': 0.0393, 'R': 0.0553, 'S': 0.0656,
    'T': 0.0534, 'V': 0.0687, 'W': 0.0108, 'Y': 0.0292
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def benjamini_hochberg(p_values):
    """Apply FDR correction."""
    n = len(p_values)
    if n == 0:
        return []
    sorted_pairs = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0] * n
    prev_adj = 0
    for rank, (orig_idx, p) in enumerate(sorted_pairs, 1):
        adj_p = min(1.0, p * n / rank)
        adj_p = max(adj_p, prev_adj)
        adjusted[orig_idx] = adj_p
        prev_adj = adj_p
    return adjusted


def calculate_enrichment(observed, total, expected_freq):
    """Calculate fold enrichment and p-value."""
    expected = expected_freq * total
    if expected == 0:
        return None, None, None

    fold = observed / expected
    try:
        result = stats.binomtest(observed, total, expected_freq, alternative='greater')
        p_value = result.pvalue
    except:
        try:
            p_value = stats.binom_test(observed, total, expected_freq, alternative='greater')
        except:
            p_value = 1.0

    return fold, p_value, expected


# =============================================================================
# MOTIF ANALYSIS FUNCTIONS
# =============================================================================

def analyze_single_aa(df, positions):
    """Analyze single amino acid enrichment at each position."""
    results = []

    for pos in positions:
        col = 'AA_site' if pos == 0 else (f'AA_+{pos}' if pos > 0 else f'AA_{pos}')
        if col not in df.columns:
            continue

        aa_counts = df[col].value_counts()
        total = df[col].notna().sum()

        if total < MIN_MOTIF_COUNT:
            continue

        for aa in AMINO_ACIDS:
            count = aa_counts.get(aa, 0)
            if count < MIN_MOTIF_COUNT:
                continue

            exp_freq = BACKGROUND_AA_FREQ.get(aa, 0.05)
            fold, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold is not None:
                results.append({
                    'Position': pos, 'Motif': aa, 'Type': 'single_AA',
                    'Observed': count, 'Expected': round(exp_count, 1), 'Total': total,
                    'Obs_freq': count / total, 'Exp_freq': exp_freq,
                    'Fold': fold, 'P_value': p_val
                })

    return pd.DataFrame(results)


def analyze_consensus_motifs(df):
    """Analyze known SUMO consensus motifs."""
    results = []

    col_m2, col_m1, col_p1, col_p2 = 'AA_-2', 'AA_-1', 'AA_+1', 'AA_+2'
    if not all(c in df.columns for c in [col_m2, col_m1, col_p1, col_p2]):
        return pd.DataFrame(results)

    valid = df[df[col_m2].notna() & df[col_m1].notna() & df[col_p1].notna() & df[col_p2].notna()]
    total = len(valid)

    if total < MIN_MOTIF_COUNT:
        return pd.DataFrame(results)

    hydrophobic = set('AVILMFYWCP')
    acidic = set('DE')

    patterns = [
        ('ψKxE (forward)', {col_m1: hydrophobic, col_p2: acidic}),
        ('ExKψ (inverse)', {col_m2: acidic, col_p1: hydrophobic}),
        ('K-x-E', {col_p2: {'E'}}),
        ('K-x-D', {col_p2: {'D'}}),
        ('K-x-[acidic]', {col_p2: acidic}),
        ('[hydrophobic]-K', {col_m1: hydrophobic}),
        ('E-x-K', {col_m2: {'E'}}),
        ('D-x-K', {col_m2: {'D'}}),
        ('[acidic]-x-K', {col_m2: acidic}),
        ('K-[hydrophobic]', {col_p1: hydrophobic}),
    ]

    for motif_name, checks in patterns:
        mask = pd.Series([True] * len(valid), index=valid.index)
        for col, valid_aas in checks.items():
            mask = mask & valid[col].isin(valid_aas)

        count = mask.sum()
        if count < MIN_MOTIF_COUNT:
            continue

        exp_freq = 1.0
        for col, valid_aas in checks.items():
            exp_freq *= sum(BACKGROUND_AA_FREQ.get(aa, 0.05) for aa in valid_aas)

        fold, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

        if fold is not None:
            results.append({
                'Position': 'consensus', 'Motif': motif_name, 'Type': 'consensus',
                'Observed': count, 'Expected': round(exp_count, 1), 'Total': total,
                'Obs_freq': count / total, 'Exp_freq': exp_freq,
                'Fold': fold, 'P_value': p_val
            })

    return pd.DataFrame(results)


def analyze_dipeptides(df, positions):
    """Analyze di-peptide enrichment."""
    results = []

    for i, pos1 in enumerate(positions[:-1]):
        pos2 = positions[i + 1]
        col1 = 'AA_site' if pos1 == 0 else (f'AA_+{pos1}' if pos1 > 0 else f'AA_{pos1}')
        col2 = 'AA_site' if pos2 == 0 else (f'AA_+{pos2}' if pos2 > 0 else f'AA_{pos2}')

        if col1 not in df.columns or col2 not in df.columns:
            continue

        valid = df[df[col1].notna() & df[col2].notna()]
        if len(valid) < MIN_MOTIF_COUNT:
            continue

        dipeptides = valid[col1].astype(str) + valid[col2].astype(str)
        counts = dipeptides.value_counts()
        total = len(dipeptides)

        for dipep, count in counts.items():
            if count < MIN_MOTIF_COUNT or len(dipep) != 2:
                continue
            if dipep[0] not in AMINO_ACIDS or dipep[1] not in AMINO_ACIDS:
                continue

            exp_freq = BACKGROUND_AA_FREQ.get(dipep[0], 0.05) * BACKGROUND_AA_FREQ.get(dipep[1], 0.05)
            fold, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold is not None and fold > 1.0:
                results.append({
                    'Position': f'{pos1},{pos2}', 'Motif': dipep, 'Type': 'dipeptide',
                    'Observed': count, 'Expected': round(exp_count, 1), 'Total': total,
                    'Obs_freq': count / total, 'Exp_freq': exp_freq,
                    'Fold': fold, 'P_value': p_val
                })

    return pd.DataFrame(results)


def run_all_analyses(df):
    """Run all motif analyses."""
    positions = list(range(-8, 9))
    all_results = []

    # Single AA
    single = analyze_single_aa(df, positions)
    if len(single) > 0:
        all_results.append(single)

    # Consensus motifs
    consensus = analyze_consensus_motifs(df)
    if len(consensus) > 0:
        all_results.append(consensus)

    # Dipeptides
    dipep = analyze_dipeptides(df, positions)
    if len(dipep) > 0:
        all_results.append(dipep)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)

    # FDR correction
    if len(combined) > 0 and 'P_value' in combined.columns:
        combined['FDR'] = benjamini_hochberg(combined['P_value'].tolist())
        combined['Significant'] = combined['FDR'] < FDR_THRESHOLD

    combined = combined.sort_values(['Type', 'Fold'], ascending=[True, False])
    return combined


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def analyze_motifs(input_file, output_file=None):
    """Perform motif discovery."""
    print(f"Analyzing motifs in {input_file}...")

    xls = pd.ExcelFile(input_file)
    sheet_name = 'Sites' if 'Sites' in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(input_file, sheet_name=sheet_name)

    print(f"Loaded {len(df)} sites")

    valid_df = df[df['pLDDT_site'].notna()].copy()
    print(f"Valid sites: {len(valid_df)}")

    score_col = 'Score (SUMO site)'
    if score_col in valid_df.columns:
        valid_df[score_col] = pd.to_numeric(valid_df[score_col], errors='coerce')

    # Define groups
    groups = {'All': valid_df}

    if score_col in valid_df.columns:
        groups['VeryHigh'] = valid_df[valid_df[score_col] >= SCORE_VERY_HIGH_THRESHOLD]
        groups['High'] = valid_df[(valid_df[score_col] >= SCORE_HIGH_THRESHOLD) & (valid_df[score_col] < SCORE_VERY_HIGH_THRESHOLD)]
        groups['Medium'] = valid_df[(valid_df[score_col] >= SCORE_MEDIUM_THRESHOLD) & (valid_df[score_col] < SCORE_HIGH_THRESHOLD)]
        groups['Low'] = valid_df[valid_df[score_col] < SCORE_MEDIUM_THRESHOLD]

    if 'Flexible' in valid_df.columns:
        groups['Flexible'] = valid_df[valid_df['Flexible'] == 'Yes']
    if 'Structured' in valid_df.columns:
        groups['Structured'] = valid_df[valid_df['Structured'] == 'Yes']

    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_motifs_simple.xlsx'

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for group_name, group_df in groups.items():
            if len(group_df) < MIN_MOTIF_COUNT:
                print(f"  Skipping {group_name} (N={len(group_df)})")
                continue

            print(f"  Analyzing {group_name} (N={len(group_df)})...")

            results = run_all_analyses(group_df)

            if len(results) > 0:
                n_sig = results['Significant'].sum() if 'Significant' in results.columns else 0
                print(f"    Found {len(results)} motifs, {n_sig} significant")
                results.to_excel(writer, sheet_name=f'Motifs_{group_name}'[:31], index=False)

    print(f"Output: {output_file}")
    return output_file


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    analyze_motifs(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
