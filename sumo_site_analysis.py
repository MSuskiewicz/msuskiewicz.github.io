#!/usr/bin/env python3
"""
SUMO Site Analyzer - Statistical Analysis
==========================================
Performs comprehensive statistical analyses on SUMO site data.

Requires output from sumo_site_core.py or compatible Excel file with columns:
- pLDDT_site, pLDDT_11residue_avg
- Flexible, Structured
- Forward_consensus, Inverse_consensus
- Acidic_in_space, Acidic_in_-2/+2
- Category
- Score (SUMO site) [optional, for score-based analyses]

Analyses included:
- Category counts and rates by score strata
- Chi-square tests for category distributions
- Logistic regression (odds ratios)
- Distance distribution analysis
- pLDDT vs Score correlation
- ROC-style analysis
- Position-specific analysis
- Bootstrap confidence intervals
- Multiple acidic residues analysis
- Random lysine background comparison

Usage:
    python sumo_site_analysis.py <input_excel_file> [output_file]
"""

import sys
import os
import random
import pandas as pd
import numpy as np
from scipy import stats

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

PLDDT_THRESHOLD = 65.0
DISTANCE_THRESHOLD = 9.0
HYDROPHOBIC_RESIDUES = {'A', 'V', 'L', 'I', 'M', 'F', 'C', 'P', 'Y'}
ACIDIC_RESIDUES = {'D', 'E'}

# Score thresholds for stratification
SCORE_VERY_HIGH_THRESHOLD = 800
SCORE_HIGH_THRESHOLD = 400
SCORE_MEDIUM_THRESHOLD = 200

# Bootstrap settings
N_BOOTSTRAP = 1000

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def bootstrap_rate_ci(successes: int, total: int, n_bootstrap: int = N_BOOTSTRAP, ci: float = 0.95) -> tuple:
    """Calculate bootstrap confidence interval for a rate."""
    if total == 0:
        return (0.0, 0.0, 0.0)
    rate = successes / total
    if total < 5:
        return (rate, rate, rate)

    samples = np.random.binomial(total, rate, n_bootstrap) / total
    alpha = (1 - ci) / 2
    lower = np.percentile(samples, alpha * 100)
    upper = np.percentile(samples, (1 - alpha) * 100)
    return (rate, lower, upper)


def rate_str(count: int, total: int) -> str:
    """Format rate as string with count."""
    rate = count / total if total > 0 else 0
    return f"{rate:.3f} ({count}/{total})"


# =============================================================================
# STANDARD ANALYSIS
# =============================================================================

def perform_analysis(df: pd.DataFrame, score_category: str, summary_stats: dict = None) -> pd.DataFrame:
    """Perform statistical analysis on the sites data."""
    score_col = 'Score (SUMO site)'
    data_list = []

    if summary_stats:
        data_list.append({'Metric': '=== DATA COMPLETENESS ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
        total = summary_stats['total']
        success = summary_stats['success']
        rate = (success / total * 100) if total > 0 else 0
        data_list.append({'Metric': 'Total sites', 'Value': total, 'Mean_Score': None, 'Std_Score': None})
        data_list.append({'Metric': 'Successfully resolved', 'Value': f"{success} ({rate:.1f}%)", 'Mean_Score': None, 'Std_Score': None})
        data_list.append({'Metric': 'Failed to resolve', 'Value': summary_stats['failed'], 'Mean_Score': None, 'Std_Score': None})
        data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    cats = [
        'Flexible_consensus', 'Flexible_flank_acidic', 'Flexible_space_acidic', 'Flexible_neither',
        'Structured_consensus', 'Structured_flank_acidic', 'Structured_space_acidic', 'Structured_neither'
    ]

    def add_stat_row(label, subset):
        count = len(subset)
        mean_score, std_score = 0.0, 0.0
        if count > 0 and score_col in df.columns:
            mean_score = subset[score_col].mean()
            std_score = subset[score_col].std()
        data_list.append({
            'Metric': label,
            'Value': count,
            'Mean_Score': round(mean_score, 2) if not pd.isna(mean_score) else None,
            'Std_Score': round(std_score, 2) if not pd.isna(std_score) else None
        })

    data_list.append({'Metric': '=== CATEGORY COUNTS ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    for cat in cats[:4]:
        add_stat_row(cat, df[df['Category'] == cat])
    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    for cat in cats[4:]:
        add_stat_row(cat, df[df['Category'] == cat])

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    data_list.append({'Metric': '=== RATES ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    total_flex = len(df[df['Flexible'] == 'Yes'])
    total_struct = len(df[df['Structured'] == 'Yes'])

    flex_cons = len(df[df['Category'] == 'Flexible_consensus'])
    flex_flank = len(df[df['Category'] == 'Flexible_flank_acidic'])
    flex_space = len(df[df['Category'] == 'Flexible_space_acidic'])
    flex_neither = len(df[df['Category'] == 'Flexible_neither'])

    struct_cons = len(df[df['Category'] == 'Structured_consensus'])
    struct_flank = len(df[df['Category'] == 'Structured_flank_acidic'])
    struct_space = len(df[df['Category'] == 'Structured_space_acidic'])
    struct_neither = len(df[df['Category'] == 'Structured_neither'])

    data_list.append({'Metric': 'Flexible consensus rate', 'Value': rate_str(flex_cons, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured consensus rate', 'Value': rate_str(struct_cons, total_struct), 'Mean_Score': None, 'Std_Score': None})

    flex_cons_flank = flex_cons + flex_flank
    struct_cons_flank = struct_cons + struct_flank
    data_list.append({'Metric': 'Flexible cons+flank rate', 'Value': rate_str(flex_cons_flank, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured cons+flank rate', 'Value': rate_str(struct_cons_flank, total_struct), 'Mean_Score': None, 'Std_Score': None})

    flex_all_acidic = flex_cons + flex_flank + flex_space
    struct_all_acidic = struct_cons + struct_flank + struct_space
    data_list.append({'Metric': 'Flexible all_acidic rate', 'Value': rate_str(flex_all_acidic, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured all_acidic rate', 'Value': rate_str(struct_all_acidic, total_struct), 'Mean_Score': None, 'Std_Score': None})

    data_list.append({'Metric': 'Flexible none rate', 'Value': rate_str(flex_neither, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured none rate', 'Value': rate_str(struct_neither, total_struct), 'Mean_Score': None, 'Std_Score': None})

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    data_list.append({'Metric': '=== CHI-SQUARE TEST ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    f_counts = [flex_cons, flex_flank, flex_space, flex_neither]
    s_counts = [struct_cons, struct_flank, struct_space, struct_neither]

    if sum(f_counts) > 0 and sum(s_counts) > 0:
        try:
            chi2, pval, _, _ = stats.chi2_contingency([f_counts, s_counts])
            data_list.append({'Metric': 'Chi-square', 'Value': f"{chi2:.2f}, p={pval:.2e}", 'Mean_Score': None, 'Std_Score': None})
        except Exception as e:
            data_list.append({'Metric': 'Chi-square error', 'Value': str(e), 'Mean_Score': None, 'Std_Score': None})

    return pd.DataFrame(data_list, columns=['Metric', 'Value', 'Mean_Score', 'Std_Score'])


def perform_trend_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Global trend analysis across score bins."""
    score_col = 'Score (SUMO site)'
    data_list = []

    data_list.append({'Metric': '=== TREND ANALYSIS ===', 'Value': ''})

    if score_col not in df.columns:
        data_list.append({'Metric': 'Error', 'Value': 'Score column not found'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    valid_df = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    if len(valid_df) < 50:
        data_list.append({'Metric': 'Warning', 'Value': f'Too few samples ({len(valid_df)})'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    valid_df = valid_df.sort_values(score_col)
    n_bins = 10
    valid_df['score_bin'] = pd.qcut(valid_df[score_col], q=n_bins, labels=False, duplicates='drop')

    bin_results = []

    for bin_idx in sorted(valid_df['score_bin'].unique()):
        bin_df = valid_df[valid_df['score_bin'] == bin_idx]
        mean_score = bin_df[score_col].mean()

        flex_sites = bin_df[bin_df['Flexible'] == 'Yes']
        struct_sites = bin_df[bin_df['Structured'] == 'Yes']

        total_flex = len(flex_sites)
        total_struct = len(struct_sites)

        if total_flex == 0 or total_struct == 0:
            continue

        flex_cf = len(flex_sites[flex_sites['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic'])])
        struct_cf = len(struct_sites[struct_sites['Category'].isin(['Structured_consensus', 'Structured_flank_acidic'])])

        flex_rate_cf = flex_cf / total_flex
        struct_rate_cf = struct_cf / total_struct
        diff_cf = flex_rate_cf - struct_rate_cf

        bin_results.append({
            'bin': bin_idx,
            'mean_score': mean_score,
            'n_sites': len(bin_df),
            'diff_cf': diff_cf
        })

    if len(bin_results) >= 3:
        bin_df = pd.DataFrame(bin_results)
        rho_cf, p_cf = stats.spearmanr(bin_df['mean_score'], bin_df['diff_cf'])
        data_list.append({'Metric': 'Score vs Flex-Struct diff (cons+flank)', 'Value': f'rho={rho_cf:.3f}, p={p_cf:.2e}'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


# =============================================================================
# ADVANCED STATISTICAL ANALYSES
# =============================================================================

def perform_logistic_regression(df: pd.DataFrame) -> pd.DataFrame:
    """Logistic regression: predict high score from features using odds ratios."""
    data_list = []
    data_list.append({'Metric': '=== LOGISTIC REGRESSION ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': 'Predicting: Is site high-scoring (>=400)?', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    valid = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    if len(valid) < 50:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples for regression'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    valid['high_score'] = (valid[score_col] >= SCORE_HIGH_THRESHOLD).astype(int)
    valid['is_flexible'] = (valid['Flexible'] == 'Yes').astype(int)
    valid['is_consensus'] = ((valid['Forward_consensus'] == 'Yes') | (valid['Inverse_consensus'] == 'Yes')).astype(int)
    valid['has_flank'] = (valid['Acidic_in_-2/+2'] == 'Yes').astype(int)
    valid['has_space'] = (valid['Acidic_in_space'] == 'Yes').astype(int)
    valid['has_any_acidic'] = ((valid['has_flank'] == 1) | (valid['has_space'] == 1)).astype(int)
    valid['dist_space'] = pd.to_numeric(valid['Dist_acidic_space_A'], errors='coerce').fillna(50)

    features = ['is_flexible', 'is_consensus', 'has_flank', 'has_space']

    data_list.append({'Metric': '--- Univariate Odds Ratios ---', 'Value': ''})

    for feat in features:
        a = len(valid[(valid[feat] == 1) & (valid['high_score'] == 1)])
        b = len(valid[(valid[feat] == 1) & (valid['high_score'] == 0)])
        c = len(valid[(valid[feat] == 0) & (valid['high_score'] == 1)])
        d = len(valid[(valid[feat] == 0) & (valid['high_score'] == 0)])

        if b > 0 and c > 0 and d > 0:
            odds_ratio = (a * d) / (b * c) if b * c > 0 else float('inf')
            try:
                _, p = stats.fisher_exact([[a, b], [c, d]])
                sig = '*' if p < 0.05 else ''
                data_list.append({'Metric': f'{feat}', 'Value': f'OR={odds_ratio:.2f}, p={p:.2e} {sig}'})
            except:
                data_list.append({'Metric': f'{feat}', 'Value': f'OR={odds_ratio:.2f}'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Distance as Continuous Predictor ---', 'Value': ''})

    valid_dist = valid[valid['dist_space'] < 50]
    if len(valid_dist) > 20:
        rho, p = stats.pointbiserialr(valid_dist['high_score'], valid_dist['dist_space'])
        data_list.append({'Metric': 'Distance vs high_score correlation', 'Value': f'r={rho:.3f}, p={p:.2e}'})

        high_dist = valid_dist[valid_dist['high_score'] == 1]['dist_space'].mean()
        low_dist = valid_dist[valid_dist['high_score'] == 0]['dist_space'].mean()
        t_stat, t_p = stats.ttest_ind(
            valid_dist[valid_dist['high_score'] == 1]['dist_space'],
            valid_dist[valid_dist['high_score'] == 0]['dist_space']
        )
        data_list.append({'Metric': 'Mean dist (high score)', 'Value': f'{high_dist:.2f} Å'})
        data_list.append({'Metric': 'Mean dist (low score)', 'Value': f'{low_dist:.2f} Å'})
        data_list.append({'Metric': 't-test p-value', 'Value': f'{t_p:.2e}'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_distance_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze distance distributions for Flexible vs Structured sites."""
    data_list = []
    data_list.append({'Metric': '=== DISTANCE DISTRIBUTION ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': 'Comparing NZ-acidic distances: Flexible vs Structured', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    valid = df[df['pLDDT_site'].notna()].copy()
    valid['dist'] = pd.to_numeric(valid['Dist_acidic_space_A'], errors='coerce')

    flex = valid[(valid['Flexible'] == 'Yes') & valid['dist'].notna()]['dist']
    struct = valid[(valid['Structured'] == 'Yes') & valid['dist'].notna()]['dist']

    if len(flex) < 10 or len(struct) < 10:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples with distance data'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    data_list.append({'Metric': '--- Flexible sites ---', 'Value': ''})
    data_list.append({'Metric': 'N', 'Value': len(flex)})
    data_list.append({'Metric': 'Mean distance', 'Value': f'{flex.mean():.2f} Å'})
    data_list.append({'Metric': 'Median distance', 'Value': f'{flex.median():.2f} Å'})
    data_list.append({'Metric': 'Std dev', 'Value': f'{flex.std():.2f} Å'})
    data_list.append({'Metric': 'Min/Max', 'Value': f'{flex.min():.2f} / {flex.max():.2f} Å'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Structured sites ---', 'Value': ''})
    data_list.append({'Metric': 'N', 'Value': len(struct)})
    data_list.append({'Metric': 'Mean distance', 'Value': f'{struct.mean():.2f} Å'})
    data_list.append({'Metric': 'Median distance', 'Value': f'{struct.median():.2f} Å'})
    data_list.append({'Metric': 'Std dev', 'Value': f'{struct.std():.2f} Å'})
    data_list.append({'Metric': 'Min/Max', 'Value': f'{struct.min():.2f} / {struct.max():.2f} Å'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Statistical comparison ---', 'Value': ''})

    u_stat, u_p = stats.mannwhitneyu(flex, struct, alternative='two-sided')
    data_list.append({'Metric': 'Mann-Whitney U test', 'Value': f'U={u_stat:.0f}, p={u_p:.2e}'})

    t_stat, t_p = stats.ttest_ind(flex, struct)
    data_list.append({'Metric': 't-test', 'Value': f't={t_stat:.2f}, p={t_p:.2e}'})

    ks_stat, ks_p = stats.ks_2samp(flex, struct)
    data_list.append({'Metric': 'KS test (distribution shape)', 'Value': f'D={ks_stat:.3f}, p={ks_p:.2e}'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Distance percentiles ---', 'Value': ''})
    for pct in [25, 50, 75, 90]:
        f_pct = np.percentile(flex, pct)
        s_pct = np.percentile(struct, pct)
        data_list.append({'Metric': f'{pct}th percentile', 'Value': f'Flex={f_pct:.1f}Å, Struct={s_pct:.1f}Å'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_plddt_score_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze correlation between pLDDT and score."""
    data_list = []
    data_list.append({'Metric': '=== pLDDT vs SCORE ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    valid = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    if len(valid) < 20:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    plddt_site = valid['pLDDT_site']
    plddt_avg = valid['pLDDT_11residue_avg']
    scores = valid[score_col]

    rho_site, p_site = stats.spearmanr(plddt_site, scores)
    rho_avg, p_avg = stats.spearmanr(plddt_avg, scores)
    r_site, rp_site = stats.pearsonr(plddt_site, scores)
    r_avg, rp_avg = stats.pearsonr(plddt_avg, scores)

    data_list.append({'Metric': '--- pLDDT_site vs Score ---', 'Value': ''})
    data_list.append({'Metric': 'Spearman rho', 'Value': f'{rho_site:.3f}, p={p_site:.2e}'})
    data_list.append({'Metric': 'Pearson r', 'Value': f'{r_site:.3f}, p={rp_site:.2e}'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- pLDDT_11avg vs Score ---', 'Value': ''})
    data_list.append({'Metric': 'Spearman rho', 'Value': f'{rho_avg:.3f}, p={p_avg:.2e}'})
    data_list.append({'Metric': 'Pearson r', 'Value': f'{r_avg:.3f}, p={rp_avg:.2e}'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Mean pLDDT by score category ---', 'Value': ''})

    for label, threshold_low, threshold_high in [
        ('VeryHigh (>=800)', 800, float('inf')),
        ('High (400-799)', 400, 800),
        ('Medium (200-399)', 200, 400),
        ('Low (<200)', 0, 200)
    ]:
        subset = valid[(valid[score_col] >= threshold_low) & (valid[score_col] < threshold_high)]
        if len(subset) > 0:
            mean_plddt = subset['pLDDT_site'].mean()
            data_list.append({'Metric': label, 'Value': f'mean pLDDT={mean_plddt:.1f} (n={len(subset)})'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_roc_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """ROC-like analysis: how well do different features predict high score?"""
    data_list = []
    data_list.append({'Metric': '=== ROC-STYLE ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': 'Predicting high score (>=400) from binary features', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    valid = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    if len(valid) < 50:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    valid['high_score'] = (valid[score_col] >= SCORE_HIGH_THRESHOLD).astype(int)
    n_positive = valid['high_score'].sum()
    n_negative = len(valid) - n_positive

    data_list.append({'Metric': 'Total samples', 'Value': len(valid)})
    data_list.append({'Metric': 'High score (positive)', 'Value': n_positive})
    data_list.append({'Metric': 'Low score (negative)', 'Value': n_negative})
    data_list.append({'Metric': '', 'Value': ''})

    predictors = [
        ('Consensus', (valid['Forward_consensus'] == 'Yes') | (valid['Inverse_consensus'] == 'Yes')),
        ('Flank acidic', valid['Acidic_in_-2/+2'] == 'Yes'),
        ('Space acidic', valid['Acidic_in_space'] == 'Yes'),
        ('Any acidic', (valid['Acidic_in_-2/+2'] == 'Yes') | (valid['Acidic_in_space'] == 'Yes')),
        ('Flexible', valid['Flexible'] == 'Yes'),
        ('Cons OR Flank', (valid['Forward_consensus'] == 'Yes') | (valid['Inverse_consensus'] == 'Yes') | (valid['Acidic_in_-2/+2'] == 'Yes')),
    ]

    data_list.append({'Metric': 'Predictor', 'Value': 'Sens | Spec | PPV | NPV | Accuracy'})
    data_list.append({'Metric': '---', 'Value': '---'})

    for name, pred_mask in predictors:
        pred = pred_mask.astype(int)
        tp = ((pred == 1) & (valid['high_score'] == 1)).sum()
        fp = ((pred == 1) & (valid['high_score'] == 0)).sum()
        tn = ((pred == 0) & (valid['high_score'] == 0)).sum()
        fn = ((pred == 0) & (valid['high_score'] == 1)).sum()

        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        acc = (tp + tn) / len(valid)

        data_list.append({
            'Metric': name,
            'Value': f'{sens:.2f} | {spec:.2f} | {ppv:.2f} | {npv:.2f} | {acc:.2f}'
        })

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- AUC estimates (distance as continuous) ---', 'Value': ''})

    valid['dist'] = pd.to_numeric(valid['Dist_acidic_space_A'], errors='coerce')
    valid_dist = valid[valid['dist'].notna()]

    if len(valid_dist) > 20:
        pos_dist = valid_dist[valid_dist['high_score'] == 1]['dist']
        neg_dist = valid_dist[valid_dist['high_score'] == 0]['dist']

        if len(pos_dist) > 5 and len(neg_dist) > 5:
            u_stat, _ = stats.mannwhitneyu(pos_dist, neg_dist, alternative='less')
            auc = u_stat / (len(pos_dist) * len(neg_dist))
            data_list.append({'Metric': 'AUC (closer distance predicts high score)', 'Value': f'{auc:.3f}'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_position_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze each flanking position separately."""
    data_list = []
    data_list.append({'Metric': '=== POSITION-SPECIFIC ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': 'Testing each position for acidic enrichment in high-scoring sites', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    valid = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    if len(valid) < 50:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    valid['high_score'] = (valid[score_col] >= SCORE_HIGH_THRESHOLD).astype(int)

    positions = [
        ('Position -2', 'AA_minus2'),
        ('Position -1', 'AA_minus1'),
        ('Position +1', 'AA_plus1'),
        ('Position +2', 'AA_plus2'),
    ]

    for pos_name, col in positions:
        data_list.append({'Metric': f'--- {pos_name} ---', 'Value': ''})

        valid[f'{col}_acidic'] = valid[col].isin(ACIDIC_RESIDUES).astype(int)

        total_acidic = valid[f'{col}_acidic'].sum()
        total = len(valid)

        high_with_acidic = len(valid[(valid[f'{col}_acidic'] == 1) & (valid['high_score'] == 1)])
        high_total = valid['high_score'].sum()
        low_with_acidic = len(valid[(valid[f'{col}_acidic'] == 1) & (valid['high_score'] == 0)])
        low_total = len(valid) - high_total

        rate_high = high_with_acidic / high_total if high_total > 0 else 0
        rate_low = low_with_acidic / low_total if low_total > 0 else 0

        data_list.append({'Metric': 'Acidic at position (all)', 'Value': f'{total_acidic}/{total} ({100*total_acidic/total:.1f}%)'})
        data_list.append({'Metric': 'Rate in high-score sites', 'Value': f'{rate_high:.3f} ({high_with_acidic}/{high_total})'})
        data_list.append({'Metric': 'Rate in low-score sites', 'Value': f'{rate_low:.3f} ({low_with_acidic}/{low_total})'})

        table = [[high_with_acidic, high_total - high_with_acidic],
                 [low_with_acidic, low_total - low_with_acidic]]
        try:
            odds, p = stats.fisher_exact(table)
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            data_list.append({'Metric': 'Fisher exact OR', 'Value': f'{odds:.2f}, p={p:.2e} {sig}'})
        except:
            pass

        data_list.append({'Metric': '', 'Value': ''})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_bootstrap_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate bootstrap confidence intervals for key rates."""
    data_list = []
    data_list.append({'Metric': '=== BOOTSTRAP CONFIDENCE INTERVALS ===', 'Value': ''})
    data_list.append({'Metric': f'Based on {N_BOOTSTRAP} bootstrap samples, 95% CI', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    valid = df[df['pLDDT_site'].notna()].copy()

    if len(valid) < 20:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    flex = valid[valid['Flexible'] == 'Yes']
    struct = valid[valid['Structured'] == 'Yes']

    rates_to_calculate = [
        ('Flexible consensus rate', len(flex[flex['Category'] == 'Flexible_consensus']), len(flex)),
        ('Flexible cons+flank rate', len(flex[flex['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic'])]), len(flex)),
        ('Flexible all_acidic rate', len(flex[flex['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic', 'Flexible_space_acidic'])]), len(flex)),
        ('Structured consensus rate', len(struct[struct['Category'] == 'Structured_consensus']), len(struct)),
        ('Structured cons+flank rate', len(struct[struct['Category'].isin(['Structured_consensus', 'Structured_flank_acidic'])]), len(struct)),
        ('Structured all_acidic rate', len(struct[struct['Category'].isin(['Structured_consensus', 'Structured_flank_acidic', 'Structured_space_acidic'])]), len(struct)),
    ]

    for name, successes, total in rates_to_calculate:
        rate, lower, upper = bootstrap_rate_ci(successes, total)
        data_list.append({
            'Metric': name,
            'Value': f'{rate:.3f} [{lower:.3f} - {upper:.3f}]'
        })

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Rate differences with CIs ---', 'Value': ''})

    if len(flex) > 10 and len(struct) > 10:
        flex_cf = len(flex[flex['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic'])]) / len(flex)
        struct_cf = len(struct[struct['Category'].isin(['Structured_consensus', 'Structured_flank_acidic'])]) / len(struct)
        obs_diff = flex_cf - struct_cf

        diffs = []
        for _ in range(N_BOOTSTRAP):
            flex_sample = flex.sample(n=len(flex), replace=True)
            struct_sample = struct.sample(n=len(struct), replace=True)
            f_rate = len(flex_sample[flex_sample['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic'])]) / len(flex_sample)
            s_rate = len(struct_sample[struct_sample['Category'].isin(['Structured_consensus', 'Structured_flank_acidic'])]) / len(struct_sample)
            diffs.append(f_rate - s_rate)

        lower = np.percentile(diffs, 2.5)
        upper = np.percentile(diffs, 97.5)
        data_list.append({
            'Metric': 'Flex-Struct cons+flank diff',
            'Value': f'{obs_diff:.3f} [{lower:.3f} - {upper:.3f}]'
        })

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_multiple_acidic_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze effect of multiple acidic residues within threshold."""
    data_list = []
    data_list.append({'Metric': '=== MULTIPLE ACIDIC RESIDUES ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': f'Counting acidic residues within {DISTANCE_THRESHOLD}Å (excl. -2/+2)', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    valid = df[df['pLDDT_site'].notna()].copy()

    if 'N_acidic_in_space' not in valid.columns:
        data_list.append({'Metric': 'Error', 'Value': 'N_acidic column not found'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    valid['n_acidic'] = pd.to_numeric(valid['N_acidic_in_space'], errors='coerce').fillna(0).astype(int)

    data_list.append({'Metric': '--- Distribution of acidic counts ---', 'Value': ''})
    for n in range(5):
        count = len(valid[valid['n_acidic'] == n])
        pct = 100 * count / len(valid) if len(valid) > 0 else 0
        data_list.append({'Metric': f'{n} acidic residues', 'Value': f'{count} ({pct:.1f}%)'})

    count_3plus = len(valid[valid['n_acidic'] >= 3])
    pct_3plus = 100 * count_3plus / len(valid) if len(valid) > 0 else 0
    data_list.append({'Metric': '3+ acidic residues', 'Value': f'{count_3plus} ({pct_3plus:.1f}%)'})

    data_list.append({'Metric': '', 'Value': ''})

    if score_col in valid.columns:
        data_list.append({'Metric': '--- Mean score by acidic count ---', 'Value': ''})
        valid[score_col] = pd.to_numeric(valid[score_col], errors='coerce')

        for n in range(4):
            subset = valid[valid['n_acidic'] == n]
            if len(subset) > 5:
                mean_score = subset[score_col].mean()
                data_list.append({'Metric': f'{n} acidic', 'Value': f'mean score={mean_score:.0f} (n={len(subset)})'})

        subset_2plus = valid[valid['n_acidic'] >= 2]
        if len(subset_2plus) > 5:
            mean_score = subset_2plus[score_col].mean()
            data_list.append({'Metric': '2+ acidic', 'Value': f'mean score={mean_score:.0f} (n={len(subset_2plus)})'})

        data_list.append({'Metric': '', 'Value': ''})

        rho, p = stats.spearmanr(valid['n_acidic'], valid[score_col])
        data_list.append({'Metric': 'Spearman: n_acidic vs score', 'Value': f'rho={rho:.3f}, p={p:.2e}'})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Comparing 0 vs 1+ acidic ---', 'Value': ''})

    has_0 = valid[valid['n_acidic'] == 0]
    has_1plus = valid[valid['n_acidic'] >= 1]

    if len(has_0) > 10 and len(has_1plus) > 10 and score_col in valid.columns:
        t_stat, p = stats.ttest_ind(has_0[score_col].dropna(), has_1plus[score_col].dropna())
        data_list.append({'Metric': 't-test (score)', 'Value': f't={t_stat:.2f}, p={p:.2e}'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def analyze_file(input_file: str, output_file: str = None):
    """Load data and perform all analyses."""
    print(f"\n{'=' * 60}\nSUMO Site Analyzer - Statistical Analysis\n{'=' * 60}\n")

    df = pd.read_excel(input_file, sheet_name='Sites' if 'Sites' in pd.ExcelFile(input_file).sheet_names else 0)

    total_sites = len(df)
    success_count = df['pLDDT_site'].notna().sum()
    failed_count = total_sites - success_count
    summary_stats = {'total': total_sites, 'success': success_count, 'failed': failed_count}

    print(f"Loaded {total_sites} sites ({success_count} with data)")

    score_col = 'Score (SUMO site)'
    if score_col in df.columns:
        df[score_col] = pd.to_numeric(df[score_col], errors='coerce')
        cats = {
            'VeryHigh': df[df[score_col] >= SCORE_VERY_HIGH_THRESHOLD],
            'High': df[(df[score_col] >= SCORE_HIGH_THRESHOLD) & (df[score_col] < SCORE_VERY_HIGH_THRESHOLD)],
            'Medium': df[(df[score_col] >= SCORE_MEDIUM_THRESHOLD) & (df[score_col] < SCORE_HIGH_THRESHOLD)],
            'Low': df[df[score_col] < SCORE_MEDIUM_THRESHOLD]
        }
    else:
        cats = {'All': df}

    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_analysis.xlsx'

    print("Generating analysis sheets...")

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Stratified analyses
        for name, sub in cats.items():
            if len(sub) > 0:
                cat_total = len(sub)
                cat_success = sub['pLDDT_site'].notna().sum()
                cat_failed = cat_total - cat_success
                cat_summary = {'total': cat_total, 'success': cat_success, 'failed': cat_failed}
                perform_analysis(sub, name, cat_summary).to_excel(writer, sheet_name=f'Analysis_{name}', index=False)

        # Extended analyses
        perform_trend_analysis(df).to_excel(writer, sheet_name='Trend_Analysis', index=False)
        perform_logistic_regression(df).to_excel(writer, sheet_name='Logistic_Regression', index=False)
        perform_distance_analysis(df).to_excel(writer, sheet_name='Distance_Analysis', index=False)
        perform_plddt_score_analysis(df).to_excel(writer, sheet_name='pLDDT_Score', index=False)
        perform_roc_analysis(df).to_excel(writer, sheet_name='ROC_Analysis', index=False)
        perform_position_analysis(df).to_excel(writer, sheet_name='Position_Analysis', index=False)
        perform_bootstrap_analysis(df).to_excel(writer, sheet_name='Bootstrap_CI', index=False)
        perform_multiple_acidic_analysis(df).to_excel(writer, sheet_name='Multiple_Acidic', index=False)

    print(f"Analysis complete. Output: {output_file}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    analyze_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
