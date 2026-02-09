#!/usr/bin/env python3
"""
SUMO Site Analyzer - Statistical Analysis (Simplified Version)
==============================================================
Simple version without threading or session pooling.
Analyzes output from sumo_site_core.py with statistical comparisons.

Usage:
    python sumo_site_analysis_simple.py <core_output_excel> [output_file]
"""

import sys
import os
import pandas as pd
import numpy as np
from scipy import stats

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

PLDDT_THRESHOLD = 65.0
DISTANCE_THRESHOLD_MIN = 4.9
DISTANCE_THRESHOLD_MAX = 8.0
SCORE_VERY_HIGH_THRESHOLD = 600
SCORE_HIGH_THRESHOLD = 400
SCORE_MEDIUM_THRESHOLD = 200

HYDROPHOBIC_RESIDUES = {'A', 'V', 'L', 'I', 'M', 'F', 'C', 'P', 'Y'}
ACIDIC_RESIDUES = {'D', 'E'}

# =============================================================================
# STATISTICAL FUNCTIONS
# =============================================================================

def count_categories(df):
    """Count sites in each 5-category hierarchical scheme."""
    suffixes = ['consensus', 'exposed_acidic_pm2', 'exposed_acidic', 'buried_acidic', 'no_acidic']
    counts = {}

    for prefix in ['Flexible', 'Structured']:
        if prefix not in df.columns:
            for suffix in suffixes:
                counts[f'{prefix}_{suffix}'] = 0
            continue

        subset = df[df[prefix] == 'Yes']
        for suffix in suffixes:
            counts[f'{prefix}_{suffix}'] = 0

        for _, row in subset.iterrows():
            # Consensus
            is_fwd = row.get('Forward_consensus') == 'Yes'
            is_inv = row.get('Inverse_consensus') == 'Yes'
            if is_fwd or is_inv:
                counts[f'{prefix}_consensus'] += 1
                continue

            # Exposed acidic in +/-2
            if row.get('Exposed_acidic_in_pm2') == 'Yes':
                counts[f'{prefix}_exposed_acidic_pm2'] += 1
                continue

            # Exposed acidic
            if row.get('Exposed_acidic_within_threshold') == 'Yes':
                counts[f'{prefix}_exposed_acidic'] += 1
                continue

            # Any acidic (buried)
            if row.get('Any_acidic_within_threshold') == 'Yes':
                counts[f'{prefix}_buried_acidic'] += 1
                continue

            # No acidic
            counts[f'{prefix}_no_acidic'] += 1

    counts['Total_Flexible'] = sum(counts[f'Flexible_{s}'] for s in suffixes)
    counts['Total_Structured'] = sum(counts[f'Structured_{s}'] for s in suffixes)
    counts['Total'] = counts['Total_Flexible'] + counts['Total_Structured']

    return counts


def calculate_rates(counts):
    """Calculate rates within Flexible and Structured."""
    suffixes = ['consensus', 'exposed_acidic_pm2', 'exposed_acidic', 'buried_acidic', 'no_acidic']
    rates = {}
    for prefix in ['Flexible', 'Structured']:
        total = counts[f'Total_{prefix}']
        for suffix in suffixes:
            key = f'{prefix}_{suffix}'
            rates[f'{key}_rate'] = counts[key] / total if total > 0 else None
    return rates


def compare_flex_vs_struct(df, label=""):
    """Chi-square test comparing Flexible vs Structured."""
    results = []
    results.append({'Metric': f'=== {label} FLEXIBLE vs STRUCTURED ===', 'Value': ''})

    suffixes = ['consensus', 'exposed_acidic_pm2', 'exposed_acidic', 'buried_acidic', 'no_acidic']
    counts = count_categories(df)
    rates = calculate_rates(counts)

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'Category', 'Value': 'Flexible | Structured'})

    for suffix in suffixes:
        f_count = counts[f'Flexible_{suffix}']
        s_count = counts[f'Structured_{suffix}']
        f_rate = rates.get(f'Flexible_{suffix}_rate', 0) or 0
        s_rate = rates.get(f'Structured_{suffix}_rate', 0) or 0
        results.append({
            'Metric': suffix,
            'Value': f'{f_count} ({f_rate:.1%}) | {s_count} ({s_rate:.1%})'
        })

    results.append({'Metric': 'TOTAL', 'Value': f"{counts['Total_Flexible']} | {counts['Total_Structured']}"})

    # Chi-square test
    flex_counts = [counts[f'Flexible_{s}'] for s in suffixes]
    struct_counts = [counts[f'Structured_{s}'] for s in suffixes]

    if sum(flex_counts) > 0 and sum(struct_counts) > 0:
        results.append({'Metric': '', 'Value': ''})
        try:
            chi2, p, dof, _ = stats.chi2_contingency([flex_counts, struct_counts])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Chi-square test', 'Value': f'χ²={chi2:.2f}, df={dof}, p={p:.2e} {sig}'})
        except Exception as e:
            results.append({'Metric': 'Chi-square error', 'Value': str(e)})

        # Fisher exact for each category
        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': '--- Fisher exact tests ---', 'Value': ''})
        for i, suffix in enumerate(suffixes):
            f_yes, f_no = flex_counts[i], sum(flex_counts) - flex_counts[i]
            s_yes, s_no = struct_counts[i], sum(struct_counts) - struct_counts[i]
            try:
                odds, p = stats.fisher_exact([[f_yes, f_no], [s_yes, s_no]])
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                results.append({'Metric': suffix, 'Value': f'OR={odds:.2f}, p={p:.2e} {sig}'})
            except:
                pass

    return results


def analyze_score_predictors(df):
    """Analyze which features predict high score."""
    results = []
    results.append({'Metric': '=== PREDICTORS OF HIGH SCORE ===', 'Value': ''})
    results.append({'Metric': f'High score >= {SCORE_HIGH_THRESHOLD}', 'Value': ''})
    results.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    if score_col not in df.columns:
        results.append({'Metric': 'Error', 'Value': 'Score column not found'})
        return results

    valid = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()
    valid[score_col] = pd.to_numeric(valid[score_col], errors='coerce')
    valid = valid[valid[score_col].notna()]

    if len(valid) < 50:
        results.append({'Metric': 'Error', 'Value': f'Too few samples ({len(valid)})'})
        return results

    valid['high_score'] = (valid[score_col] >= SCORE_HIGH_THRESHOLD).astype(int)
    n_high = valid['high_score'].sum()
    n_low = len(valid) - n_high

    results.append({'Metric': 'Total valid sites', 'Value': len(valid)})
    results.append({'Metric': 'High score', 'Value': f'{n_high} ({100*n_high/len(valid):.1f}%)'})
    results.append({'Metric': 'Low score', 'Value': f'{n_low} ({100*n_low/len(valid):.1f}%)'})
    results.append({'Metric': '', 'Value': ''})

    # Test predictors
    predictors = [
        ('Flexible', valid.get('Flexible') == 'Yes'),
        ('Structured', valid.get('Structured') == 'Yes'),
        ('Any_consensus', (valid.get('Forward_consensus') == 'Yes') | (valid.get('Inverse_consensus') == 'Yes')),
        ('Any_acidic_within_threshold', valid.get('Any_acidic_within_threshold') == 'Yes'),
        ('Exposed_acidic_within_threshold', valid.get('Exposed_acidic_within_threshold') == 'Yes'),
        ('Exposed_acidic_in_pm2', valid.get('Exposed_acidic_in_pm2') == 'Yes'),
    ]

    results.append({'Metric': '--- Binary Predictors ---', 'Value': ''})
    results.append({'Metric': 'Predictor', 'Value': 'Rate_High | Rate_Low | OR | p-value'})

    for name, mask in predictors:
        if mask is None:
            continue
        pred = mask.astype(int)
        tp = ((pred == 1) & (valid['high_score'] == 1)).sum()
        fp = ((pred == 1) & (valid['high_score'] == 0)).sum()
        fn = ((pred == 0) & (valid['high_score'] == 1)).sum()
        tn = ((pred == 0) & (valid['high_score'] == 0)).sum()

        rate_high = tp / n_high if n_high > 0 else 0
        rate_low = fp / n_low if n_low > 0 else 0

        if fp > 0 and fn > 0:
            odds = (tp * tn) / (fp * fn) if fp * fn > 0 else float('inf')
            try:
                _, p = stats.fisher_exact([[tp, fp], [fn, tn]])
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                results.append({
                    'Metric': name,
                    'Value': f'{rate_high:.1%} | {rate_low:.1%} | {odds:.2f} | {p:.2e} {sig}'
                })
            except:
                results.append({'Metric': name, 'Value': f'{rate_high:.1%} | {rate_low:.1%} | N/A'})

    return results


def analyze_aa_enrichment(df, score_col='Score (SUMO site)'):
    """Analyze amino acid enrichment at each position."""
    positions = list(range(-8, 9))
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    pos_cols = [f'AA_{p:+d}' if p != 0 else 'AA_site' for p in positions]

    df = df.copy()
    if score_col in df.columns:
        df[score_col] = pd.to_numeric(df[score_col], errors='coerce')

    subsets = {
        'All': df,
        'Flexible': df[df['Flexible'] == 'Yes'] if 'Flexible' in df.columns else pd.DataFrame(),
        'Structured': df[df['Structured'] == 'Yes'] if 'Structured' in df.columns else pd.DataFrame(),
    }

    if score_col in df.columns:
        subsets['VeryHigh'] = df[df[score_col] >= SCORE_VERY_HIGH_THRESHOLD]
        subsets['High'] = df[(df[score_col] >= SCORE_HIGH_THRESHOLD) & (df[score_col] < SCORE_VERY_HIGH_THRESHOLD)]
        subsets['Medium'] = df[(df[score_col] >= SCORE_MEDIUM_THRESHOLD) & (df[score_col] < SCORE_HIGH_THRESHOLD)]
        subsets['Low'] = df[df[score_col] < SCORE_MEDIUM_THRESHOLD]

    all_rows = []

    for subset_name, subset_df in subsets.items():
        if len(subset_df) == 0:
            continue

        all_rows.append({'AA': f'=== {subset_name.upper()} (N={len(subset_df)}) ===', **{col: '' for col in pos_cols}})

        for aa in amino_acids:
            row = {'AA': aa}
            for i, pos in enumerate(positions):
                col_name = pos_cols[i]
                if col_name in subset_df.columns:
                    count = (subset_df[col_name] == aa).sum()
                    total = subset_df[col_name].notna().sum()
                    rate = count / total if total > 0 else 0
                    row[col_name] = f'{count} ({rate:.1%})'
                else:
                    row[col_name] = 'N/A'
            all_rows.append(row)

        all_rows.append({'AA': '', **{col: '' for col in pos_cols}})

    return pd.DataFrame(all_rows)


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def analyze_file(input_file, output_file=None):
    """Perform all analyses."""
    print(f"Analyzing {input_file}...")

    xls = pd.ExcelFile(input_file)
    sheet_name = 'Sites' if 'Sites' in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(input_file, sheet_name=sheet_name)

    print(f"Loaded {len(df)} sites")

    valid_df = df[df['pLDDT_site'].notna()].copy()
    print(f"Valid sites: {len(valid_df)}")

    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_analysis_simple.xlsx'

    score_col = 'Score (SUMO site)'
    if score_col in valid_df.columns:
        valid_df[score_col] = pd.to_numeric(valid_df[score_col], errors='coerce')

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Global analysis
        global_results = []
        global_results.append({'Metric': '=== GLOBAL ANALYSIS ===', 'Value': ''})
        global_results.append({'Metric': '', 'Value': ''})

        counts = count_categories(valid_df)
        rates = calculate_rates(counts)

        suffixes = ['consensus', 'exposed_acidic_pm2', 'exposed_acidic', 'buried_acidic', 'no_acidic']

        global_results.append({'Metric': '--- Category Counts ---', 'Value': ''})
        for suffix in suffixes:
            cat = f'Flexible_{suffix}'
            rate = rates.get(f'{cat}_rate', 0) or 0
            global_results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%})'})

        global_results.append({'Metric': '', 'Value': ''})
        for suffix in suffixes:
            cat = f'Structured_{suffix}'
            rate = rates.get(f'{cat}_rate', 0) or 0
            global_results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%})'})

        global_results.append({'Metric': '', 'Value': ''})
        global_results.append({'Metric': 'Total Flexible', 'Value': counts['Total_Flexible']})
        global_results.append({'Metric': 'Total Structured', 'Value': counts['Total_Structured']})

        global_results.append({'Metric': '', 'Value': ''})
        global_results.extend(compare_flex_vs_struct(valid_df, "ALL SITES"))

        pd.DataFrame(global_results).to_excel(writer, sheet_name='Global_Analysis', index=False)

        # Score-stratified analysis
        if score_col in valid_df.columns:
            strata = {
                'VeryHigh': valid_df[valid_df[score_col] >= SCORE_VERY_HIGH_THRESHOLD],
                'High': valid_df[(valid_df[score_col] >= SCORE_HIGH_THRESHOLD) & (valid_df[score_col] < SCORE_VERY_HIGH_THRESHOLD)],
                'Medium': valid_df[(valid_df[score_col] >= SCORE_MEDIUM_THRESHOLD) & (valid_df[score_col] < SCORE_HIGH_THRESHOLD)],
                'Low': valid_df[valid_df[score_col] < SCORE_MEDIUM_THRESHOLD]
            }

            for name, subset in strata.items():
                if len(subset) == 0:
                    continue

                results = []
                results.append({'Metric': f'=== {name.upper()} SCORE ===', 'Value': ''})
                results.append({'Metric': 'N sites', 'Value': len(subset)})
                results.append({'Metric': '', 'Value': ''})

                s_counts = count_categories(subset)
                s_rates = calculate_rates(s_counts)

                results.append({'Metric': '--- Category Counts ---', 'Value': ''})
                for suffix in suffixes:
                    cat = f'Flexible_{suffix}'
                    rate = s_rates.get(f'{cat}_rate', 0) or 0
                    results.append({'Metric': cat, 'Value': f'{s_counts[cat]} ({rate:.1%})'})

                results.append({'Metric': '', 'Value': ''})
                for suffix in suffixes:
                    cat = f'Structured_{suffix}'
                    rate = s_rates.get(f'{cat}_rate', 0) or 0
                    results.append({'Metric': cat, 'Value': f'{s_counts[cat]} ({rate:.1%})'})

                results.append({'Metric': '', 'Value': ''})
                results.extend(compare_flex_vs_struct(subset, name.upper()))

                pd.DataFrame(results).to_excel(writer, sheet_name=f'Score_{name}', index=False)

        # Score predictors
        predictors = analyze_score_predictors(valid_df)
        pd.DataFrame(predictors).to_excel(writer, sheet_name='Score_Predictors', index=False)

        # AA enrichment
        aa_enrich = analyze_aa_enrichment(valid_df)
        aa_enrich.to_excel(writer, sheet_name='AA_Enrichment', index=False)

    print(f"Analysis complete: {output_file}")
    return output_file


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    analyze_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
