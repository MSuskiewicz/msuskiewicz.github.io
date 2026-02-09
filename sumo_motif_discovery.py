#!/usr/bin/env python3
"""
SUMO Site Analyzer - De Novo Motif Discovery
=============================================
Analyzes output from sumo_site_core.py to discover enriched sequence motifs
surrounding SUMO modification sites.

Analyses performed:
1. Position-specific amino acid enrichment (single residues)
2. Di-peptide (2-mer) motif enrichment at each position
3. Tri-peptide (3-mer) motif enrichment at each position
4. Gapped motifs (e.g., X-x-X patterns)

All analyses are stratified by score groups (VeryHigh, High, Medium, Low)
and compared against background frequencies.

Usage:
    python sumo_motif_discovery.py <core_output_excel> [output_file]
"""

import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
from collections import Counter
from itertools import product
import warnings

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

# Score thresholds (must match core script)
SCORE_VERY_HIGH_THRESHOLD = 600
SCORE_HIGH_THRESHOLD = 400
SCORE_MEDIUM_THRESHOLD = 200

# Motif discovery parameters
MIN_MOTIF_COUNT = 5          # Minimum occurrences to consider a motif
P_VALUE_THRESHOLD = 0.05     # Significance threshold (before correction)
FDR_THRESHOLD = 0.1          # False discovery rate threshold

# Amino acid categories for pattern analysis
AA_CATEGORIES = {
    'hydrophobic': set('AVILMFYWCP'),
    'acidic': set('DE'),
    'basic': set('KRH'),
    'polar': set('STNQ'),
    'small': set('AGST'),
    'aromatic': set('FYW'),
    'aliphatic': set('AVILM'),
}

# Standard amino acids
AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

# Background amino acid frequencies (from UniProt/SwissProt)
# These are approximate proteome-wide frequencies
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

def get_sequence_context(row, positions):
    """Extract sequence context from a row of the DataFrame."""
    seq = []
    for pos in positions:
        if pos == 0:
            col = 'AA_site'
        elif pos > 0:
            col = f'AA_+{pos}'
        else:
            col = f'AA_{pos}'

        aa = row.get(col, None)
        if pd.isna(aa) or aa not in AMINO_ACIDS:
            seq.append(None)
        else:
            seq.append(aa)
    return seq


def calculate_enrichment(observed_count, total_observed, expected_freq, total_expected=None):
    """Calculate fold enrichment and p-value using binomial test."""
    if total_expected is None:
        total_expected = total_observed

    expected_count = expected_freq * total_observed

    if expected_count == 0:
        return None, None, None

    fold_enrichment = observed_count / expected_count if expected_count > 0 else float('inf')

    # Binomial test
    try:
        p_value = stats.binom_test(observed_count, total_observed, expected_freq,
                                    alternative='greater')
    except:
        try:
            # scipy >= 1.7 uses binomtest
            result = stats.binomtest(observed_count, total_observed, expected_freq,
                                     alternative='greater')
            p_value = result.pvalue
        except:
            p_value = 1.0

    return fold_enrichment, p_value, expected_count


def benjamini_hochberg(p_values):
    """Apply Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return []

    # Sort p-values and track original indices
    sorted_pairs = sorted(enumerate(p_values), key=lambda x: x[1])

    # Calculate adjusted p-values
    adjusted = [0] * n
    prev_adj = 0

    for rank, (orig_idx, p) in enumerate(sorted_pairs, 1):
        adj_p = min(1.0, p * n / rank)
        adj_p = max(adj_p, prev_adj)  # Ensure monotonicity
        adjusted[orig_idx] = adj_p
        prev_adj = adj_p

    return adjusted


def categorize_aa(aa):
    """Categorize amino acid into functional groups."""
    categories = []
    for cat_name, cat_aas in AA_CATEGORIES.items():
        if aa in cat_aas:
            categories.append(cat_name)
    return categories


# =============================================================================
# MOTIF ANALYSIS FUNCTIONS
# =============================================================================

def analyze_single_aa_enrichment(df, positions, background_freq=None):
    """Analyze single amino acid enrichment at each position.

    Returns DataFrame with enrichment statistics for each AA at each position.
    """
    if background_freq is None:
        background_freq = BACKGROUND_AA_FREQ

    results = []

    for pos in positions:
        if pos == 0:
            col = 'AA_site'
        elif pos > 0:
            col = f'AA_+{pos}'
        else:
            col = f'AA_{pos}'

        if col not in df.columns:
            continue

        # Count amino acids at this position
        aa_counts = df[col].value_counts()
        total = df[col].notna().sum()

        if total < MIN_MOTIF_COUNT:
            continue

        for aa in AMINO_ACIDS:
            count = aa_counts.get(aa, 0)
            if count < MIN_MOTIF_COUNT:
                continue

            obs_freq = count / total
            exp_freq = background_freq.get(aa, 0.05)

            fold_enrich, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold_enrich is not None:
                results.append({
                    'Position': pos,
                    'Motif': aa,
                    'Type': 'single_AA',
                    'Observed_count': count,
                    'Expected_count': round(exp_count, 1),
                    'Total': total,
                    'Observed_freq': obs_freq,
                    'Expected_freq': exp_freq,
                    'Fold_enrichment': fold_enrich,
                    'P_value': p_val,
                    'Categories': ','.join(categorize_aa(aa))
                })

    return pd.DataFrame(results)


def analyze_dipeptide_enrichment(df, positions):
    """Analyze di-peptide (2-mer) motif enrichment.

    Looks at consecutive amino acid pairs at each position.
    """
    results = []

    # Calculate background dipeptide frequencies
    bg_dipeptide = {}
    for aa1 in AMINO_ACIDS:
        for aa2 in AMINO_ACIDS:
            bg_dipeptide[aa1 + aa2] = BACKGROUND_AA_FREQ.get(aa1, 0.05) * BACKGROUND_AA_FREQ.get(aa2, 0.05)

    for i, pos1 in enumerate(positions[:-1]):
        pos2 = positions[i + 1]

        # Get column names
        if pos1 == 0:
            col1 = 'AA_site'
        elif pos1 > 0:
            col1 = f'AA_+{pos1}'
        else:
            col1 = f'AA_{pos1}'

        if pos2 == 0:
            col2 = 'AA_site'
        elif pos2 > 0:
            col2 = f'AA_+{pos2}'
        else:
            col2 = f'AA_{pos2}'

        if col1 not in df.columns or col2 not in df.columns:
            continue

        # Create dipeptide column
        valid_mask = df[col1].notna() & df[col2].notna()
        valid_df = df[valid_mask]

        if len(valid_df) < MIN_MOTIF_COUNT:
            continue

        dipeptides = valid_df[col1].astype(str) + valid_df[col2].astype(str)
        dipeptide_counts = dipeptides.value_counts()
        total = len(dipeptides)

        for dipep, count in dipeptide_counts.items():
            if count < MIN_MOTIF_COUNT:
                continue
            if len(dipep) != 2 or dipep[0] not in AMINO_ACIDS or dipep[1] not in AMINO_ACIDS:
                continue

            exp_freq = bg_dipeptide.get(dipep, 0.0025)
            fold_enrich, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold_enrich is not None and fold_enrich > 1.0:
                results.append({
                    'Position': f'{pos1},{pos2}',
                    'Motif': dipep,
                    'Type': 'dipeptide',
                    'Observed_count': count,
                    'Expected_count': round(exp_count, 1),
                    'Total': total,
                    'Observed_freq': count / total,
                    'Expected_freq': exp_freq,
                    'Fold_enrichment': fold_enrich,
                    'P_value': p_val,
                    'Categories': ''
                })

    return pd.DataFrame(results)


def analyze_tripeptide_enrichment(df, positions):
    """Analyze tri-peptide (3-mer) motif enrichment.

    Looks at consecutive amino acid triplets at each position.
    """
    results = []

    # Calculate background tripeptide frequencies
    bg_tripeptide = {}
    for aa1, aa2, aa3 in product(AMINO_ACIDS, repeat=3):
        tripep = aa1 + aa2 + aa3
        bg_tripeptide[tripep] = (BACKGROUND_AA_FREQ.get(aa1, 0.05) *
                                  BACKGROUND_AA_FREQ.get(aa2, 0.05) *
                                  BACKGROUND_AA_FREQ.get(aa3, 0.05))

    for i, pos1 in enumerate(positions[:-2]):
        pos2 = positions[i + 1]
        pos3 = positions[i + 2]

        # Get column names
        cols = []
        for pos in [pos1, pos2, pos3]:
            if pos == 0:
                cols.append('AA_site')
            elif pos > 0:
                cols.append(f'AA_+{pos}')
            else:
                cols.append(f'AA_{pos}')

        if not all(c in df.columns for c in cols):
            continue

        # Create tripeptide column
        valid_mask = df[cols[0]].notna() & df[cols[1]].notna() & df[cols[2]].notna()
        valid_df = df[valid_mask]

        if len(valid_df) < MIN_MOTIF_COUNT:
            continue

        tripeptides = (valid_df[cols[0]].astype(str) +
                       valid_df[cols[1]].astype(str) +
                       valid_df[cols[2]].astype(str))
        tripeptide_counts = tripeptides.value_counts()
        total = len(tripeptides)

        for tripep, count in tripeptide_counts.items():
            if count < MIN_MOTIF_COUNT:
                continue
            if len(tripep) != 3 or not all(aa in AMINO_ACIDS for aa in tripep):
                continue

            exp_freq = bg_tripeptide.get(tripep, 0.000125)
            fold_enrich, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold_enrich is not None and fold_enrich > 1.0:
                results.append({
                    'Position': f'{pos1},{pos2},{pos3}',
                    'Motif': tripep,
                    'Type': 'tripeptide',
                    'Observed_count': count,
                    'Expected_count': round(exp_count, 1),
                    'Total': total,
                    'Observed_freq': count / total,
                    'Expected_freq': exp_freq,
                    'Fold_enrichment': fold_enrich,
                    'P_value': p_val,
                    'Categories': ''
                })

    return pd.DataFrame(results)


def analyze_gapped_motifs(df, positions):
    """Analyze gapped motifs (patterns with variable positions).

    Looks for patterns like X-x-Y where x is any amino acid.
    Focuses on positions around the lysine site.
    """
    results = []

    # Define interesting gap patterns to check
    # Format: (pos1, pos2) tuples for non-adjacent positions
    gap_patterns = [
        (-2, 0),   # X-x-K
        (-2, 2),   # X-x-K-x-Y (flanking the site)
        (-1, 1),   # X-K-Y (immediate neighbors)
        (-3, -1),  # X-x-Y-K
        (1, 3),    # K-X-x-Y
        (-4, -2),
        (2, 4),
    ]

    for pos1, pos2 in gap_patterns:
        if pos1 not in positions or pos2 not in positions:
            continue

        # Get column names
        if pos1 == 0:
            col1 = 'AA_site'
        elif pos1 > 0:
            col1 = f'AA_+{pos1}'
        else:
            col1 = f'AA_{pos1}'

        if pos2 == 0:
            col2 = 'AA_site'
        elif pos2 > 0:
            col2 = f'AA_+{pos2}'
        else:
            col2 = f'AA_{pos2}'

        if col1 not in df.columns or col2 not in df.columns:
            continue

        valid_mask = df[col1].notna() & df[col2].notna()
        valid_df = df[valid_mask]

        if len(valid_df) < MIN_MOTIF_COUNT:
            continue

        # Create gapped motif
        gap_size = abs(pos2 - pos1) - 1
        gapped = (valid_df[col1].astype(str) +
                  '-' + 'x' * gap_size + '-' +
                  valid_df[col2].astype(str))

        # Also create simple version for counting
        simple_gapped = valid_df[col1].astype(str) + '_' + valid_df[col2].astype(str)
        gapped_counts = simple_gapped.value_counts()
        total = len(simple_gapped)

        for motif_key, count in gapped_counts.items():
            if count < MIN_MOTIF_COUNT:
                continue

            parts = motif_key.split('_')
            if len(parts) != 2:
                continue
            aa1, aa2 = parts
            if aa1 not in AMINO_ACIDS or aa2 not in AMINO_ACIDS:
                continue

            # Expected frequency assuming independence
            exp_freq = BACKGROUND_AA_FREQ.get(aa1, 0.05) * BACKGROUND_AA_FREQ.get(aa2, 0.05)
            fold_enrich, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold_enrich is not None and fold_enrich > 1.5:  # Higher threshold for gapped
                # Create readable motif representation
                gap_str = 'x' * gap_size
                motif_repr = f'{aa1}-{gap_str}-{aa2}'

                results.append({
                    'Position': f'{pos1}...{pos2}',
                    'Motif': motif_repr,
                    'Type': 'gapped',
                    'Observed_count': count,
                    'Expected_count': round(exp_count, 1),
                    'Total': total,
                    'Observed_freq': count / total,
                    'Expected_freq': exp_freq,
                    'Fold_enrichment': fold_enrich,
                    'P_value': p_val,
                    'Categories': ''
                })

    return pd.DataFrame(results)


def analyze_category_patterns(df, positions):
    """Analyze patterns using amino acid categories (hydrophobic, acidic, etc.)."""
    results = []

    # Check for category enrichment at each position
    for pos in positions:
        if pos == 0:
            col = 'AA_site'
        elif pos > 0:
            col = f'AA_+{pos}'
        else:
            col = f'AA_{pos}'

        if col not in df.columns:
            continue

        valid_df = df[df[col].notna()]
        total = len(valid_df)

        if total < MIN_MOTIF_COUNT:
            continue

        for cat_name, cat_aas in AA_CATEGORIES.items():
            # Count amino acids in this category
            count = valid_df[col].isin(cat_aas).sum()

            if count < MIN_MOTIF_COUNT:
                continue

            # Expected frequency based on background
            exp_freq = sum(BACKGROUND_AA_FREQ.get(aa, 0.05) for aa in cat_aas)
            fold_enrich, p_val, exp_count = calculate_enrichment(count, total, exp_freq)

            if fold_enrich is not None and fold_enrich > 1.0:
                results.append({
                    'Position': pos,
                    'Motif': f'[{cat_name}]',
                    'Type': 'category',
                    'Observed_count': count,
                    'Expected_count': round(exp_count, 1),
                    'Total': total,
                    'Observed_freq': count / total,
                    'Expected_freq': exp_freq,
                    'Fold_enrichment': fold_enrich,
                    'P_value': p_val,
                    'Categories': cat_name
                })

    return pd.DataFrame(results)


def run_all_analyses(df, label=''):
    """Run all motif analyses on a DataFrame subset."""
    positions = list(range(-8, 9))  # -8 to +8

    all_results = []

    # Single amino acid enrichment
    single_aa = analyze_single_aa_enrichment(df, positions)
    if len(single_aa) > 0:
        all_results.append(single_aa)

    # Dipeptide enrichment
    dipeptide = analyze_dipeptide_enrichment(df, positions)
    if len(dipeptide) > 0:
        all_results.append(dipeptide)

    # Tripeptide enrichment
    tripeptide = analyze_tripeptide_enrichment(df, positions)
    if len(tripeptide) > 0:
        all_results.append(tripeptide)

    # Gapped motifs
    gapped = analyze_gapped_motifs(df, positions)
    if len(gapped) > 0:
        all_results.append(gapped)

    # Category patterns
    category = analyze_category_patterns(df, positions)
    if len(category) > 0:
        all_results.append(category)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)

    # Apply FDR correction
    if len(combined) > 0 and 'P_value' in combined.columns:
        combined['FDR'] = benjamini_hochberg(combined['P_value'].tolist())
        combined['Significant'] = combined['FDR'] < FDR_THRESHOLD

    # Sort by fold enrichment
    combined = combined.sort_values(['Type', 'Fold_enrichment'], ascending=[True, False])

    return combined


def create_summary_sheet(results_dict):
    """Create a summary sheet with top motifs from each group."""
    summary = []

    summary.append({'Group': '=== TOP ENRICHED MOTIFS SUMMARY ===',
                    'Motif': '', 'Position': '', 'Fold': '', 'FDR': ''})
    summary.append({'Group': '', 'Motif': '', 'Position': '', 'Fold': '', 'FDR': ''})

    for group_name, df in results_dict.items():
        if len(df) == 0:
            continue

        # Get significant motifs
        sig_df = df[df['Significant'] == True] if 'Significant' in df.columns else df

        if len(sig_df) == 0:
            sig_df = df.head(10)  # Show top 10 even if not significant

        summary.append({'Group': f'--- {group_name} (N significant = {len(sig_df)}) ---',
                        'Motif': '', 'Position': '', 'Fold': '', 'FDR': ''})

        # Top 5 by fold enrichment for each type
        for motif_type in ['single_AA', 'dipeptide', 'tripeptide', 'gapped', 'category']:
            type_df = sig_df[sig_df['Type'] == motif_type].head(5)

            for _, row in type_df.iterrows():
                summary.append({
                    'Group': group_name,
                    'Motif': row['Motif'],
                    'Position': row['Position'],
                    'Fold': f"{row['Fold_enrichment']:.2f}",
                    'FDR': f"{row.get('FDR', row['P_value']):.2e}"
                })

        summary.append({'Group': '', 'Motif': '', 'Position': '', 'Fold': '', 'FDR': ''})

    return pd.DataFrame(summary)


def create_consensus_motif(df, positions):
    """Generate a consensus motif representation from enriched positions."""
    consensus = []

    for pos in positions:
        if pos == 0:
            consensus.append('K')  # Lysine site
            continue

        col = f'AA_+{pos}' if pos > 0 else f'AA_{pos}'
        if col not in df.columns:
            consensus.append('x')
            continue

        aa_counts = df[col].value_counts()
        total = df[col].notna().sum()

        if total == 0:
            consensus.append('x')
            continue

        # Find most enriched
        best_aa = 'x'
        best_fold = 1.0

        for aa in AMINO_ACIDS:
            count = aa_counts.get(aa, 0)
            if count < 3:
                continue

            obs_freq = count / total
            exp_freq = BACKGROUND_AA_FREQ.get(aa, 0.05)
            fold = obs_freq / exp_freq if exp_freq > 0 else 1.0

            if fold > best_fold and fold > 1.5:
                best_fold = fold
                best_aa = aa

        # Check for category enrichment
        if best_aa == 'x':
            for cat_name, cat_aas in [('hydrophobic', 'ΨVILMFYWCP'),
                                       ('acidic', 'DE'),
                                       ('basic', 'KRH')]:
                cat_count = df[col].isin(set(cat_aas)).sum()
                if cat_count / total > 0.3:
                    best_aa = f'[{cat_name[:3]}]'
                    break

        consensus.append(best_aa if best_aa != 'x' else 'x')

    return ''.join(consensus)


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def analyze_motifs(input_file: str, output_file: str = None):
    """Load core output and perform motif discovery analysis."""
    print(f"\n{'=' * 60}\nSUMO Site Analyzer - De Novo Motif Discovery\n{'=' * 60}\n")
    print(f"Parameters:")
    print(f"  Min motif count: {MIN_MOTIF_COUNT}")
    print(f"  FDR threshold: {FDR_THRESHOLD}")
    print(f"  Score thresholds: {SCORE_VERY_HIGH_THRESHOLD}/{SCORE_HIGH_THRESHOLD}/{SCORE_MEDIUM_THRESHOLD}")
    print()

    # Load data
    xls = pd.ExcelFile(input_file)
    sheet_name = 'Sites' if 'Sites' in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(input_file, sheet_name=sheet_name)

    print(f"Loaded {len(df)} sites from {input_file}")

    # Filter to valid sites
    valid_df = df[df['pLDDT_site'].notna()].copy()
    print(f"Valid sites with structure data: {len(valid_df)}")

    # Prepare score column
    score_col = 'Score (SUMO site)'
    if score_col in valid_df.columns:
        valid_df[score_col] = pd.to_numeric(valid_df[score_col], errors='coerce')

    # Define score groups
    score_groups = {
        'All': valid_df,
        'VeryHigh': valid_df[valid_df[score_col] >= SCORE_VERY_HIGH_THRESHOLD] if score_col in valid_df.columns else pd.DataFrame(),
        'High': valid_df[(valid_df[score_col] >= SCORE_HIGH_THRESHOLD) & (valid_df[score_col] < SCORE_VERY_HIGH_THRESHOLD)] if score_col in valid_df.columns else pd.DataFrame(),
        'Medium': valid_df[(valid_df[score_col] >= SCORE_MEDIUM_THRESHOLD) & (valid_df[score_col] < SCORE_HIGH_THRESHOLD)] if score_col in valid_df.columns else pd.DataFrame(),
        'Low': valid_df[valid_df[score_col] < SCORE_MEDIUM_THRESHOLD] if score_col in valid_df.columns else pd.DataFrame(),
    }

    # Also analyze by flexibility
    if 'Flexible' in valid_df.columns:
        score_groups['Flexible'] = valid_df[valid_df['Flexible'] == 'Yes']
    if 'Structured' in valid_df.columns:
        score_groups['Structured'] = valid_df[valid_df['Structured'] == 'Yes']

    # Prepare output
    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_motifs.xlsx'

    print("\nRunning motif discovery analyses...")

    results_dict = {}

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for group_name, group_df in score_groups.items():
            if len(group_df) < MIN_MOTIF_COUNT:
                print(f"  Skipping {group_name} (N={len(group_df)} < {MIN_MOTIF_COUNT})")
                continue

            print(f"  Analyzing {group_name} (N={len(group_df)})...")

            results = run_all_analyses(group_df, group_name)
            results_dict[group_name] = results

            if len(results) > 0:
                # Reorder columns for output
                col_order = ['Type', 'Position', 'Motif', 'Observed_count', 'Expected_count',
                             'Total', 'Observed_freq', 'Expected_freq', 'Fold_enrichment',
                             'P_value', 'FDR', 'Significant', 'Categories']
                col_order = [c for c in col_order if c in results.columns]
                results = results[col_order]

                # Truncate sheet name to Excel limit (31 chars)
                sheet_name = f'Motifs_{group_name}'[:31]
                results.to_excel(writer, sheet_name=sheet_name, index=False)

                # Count significant
                n_sig = results['Significant'].sum() if 'Significant' in results.columns else 0
                print(f"    Found {len(results)} motifs, {n_sig} significant (FDR<{FDR_THRESHOLD})")

        # Create summary sheet
        summary_df = create_summary_sheet(results_dict)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)

        # Create consensus motif sheet
        positions = list(range(-8, 9))
        consensus_rows = []
        consensus_rows.append({'Group': 'CONSENSUS MOTIFS', 'Consensus': ''})
        consensus_rows.append({'Group': 'Position', 'Consensus': ' '.join([f'{p:+d}' if p != 0 else ' K' for p in positions])})
        consensus_rows.append({'Group': '', 'Consensus': ''})

        for group_name, group_df in score_groups.items():
            if len(group_df) >= MIN_MOTIF_COUNT:
                consensus = create_consensus_motif(group_df, positions)
                consensus_rows.append({'Group': group_name, 'Consensus': consensus})

        pd.DataFrame(consensus_rows).to_excel(writer, sheet_name='Consensus', index=False)

    print(f"\nMotif discovery complete. Output: {output_file}")
    return output_file


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nConfigurable parameters (edit at top of script):")
        print(f"  MIN_MOTIF_COUNT = {MIN_MOTIF_COUNT}")
        print(f"  P_VALUE_THRESHOLD = {P_VALUE_THRESHOLD}")
        print(f"  FDR_THRESHOLD = {FDR_THRESHOLD}")
        print(f"  SCORE_VERY_HIGH_THRESHOLD = {SCORE_VERY_HIGH_THRESHOLD}")
        print(f"  SCORE_HIGH_THRESHOLD = {SCORE_HIGH_THRESHOLD}")
        print(f"  SCORE_MEDIUM_THRESHOLD = {SCORE_MEDIUM_THRESHOLD}")
        sys.exit(1)
    analyze_motifs(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
