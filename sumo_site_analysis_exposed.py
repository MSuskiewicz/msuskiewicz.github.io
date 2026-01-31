#!/usr/bin/env python3
"""
SUMO Site Analyzer - Statistical Analysis (Exposed Acidic Focus)
=================================================================
Analyzes output from sumo_site_core_exposed.py with statistical comparisons
focusing on surface-exposed acidic residues.

Analyses:
1. Global category counts (Flexible/Structured × exposed/noexposed)
2. Score-stratified analysis (VeryHigh, High, Medium, Low)
3. Control analysis: all lysines in analyzed AlphaFold models
4. Statistical comparisons: Flexible vs Structured, Sites vs Control
5. Predictors of high score (focusing on exposed acidics)

Usage:
    python sumo_site_analysis_exposed.py <core_output_excel> [output_file]
"""

import sys
import os
import re
import math
import requests
import pandas as pd
import numpy as np
from scipy import stats
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
import logging

# =============================================================================
# CONFIGURABLE PARAMETERS (must match sumo_site_core_exposed.py)
# =============================================================================

PLDDT_THRESHOLD = 65.0          # pLDDT below this = Flexible, above = Structured
PLDDT_WINDOW_LENGTH = 11        # Window size for pLDDT averaging (centered on site)

# Distance thresholds (separate for flexible and structured sites)
DISTANCE_THRESHOLD_FLEXIBLE = 8.0    # Angstroms for Cα-Cα distance in flexible sites
DISTANCE_THRESHOLD_STRUCTURED = 8.0  # Angstroms for NZ-CG/CD distance in structured sites

# Surface exposure parameters
EXPOSURE_NEIGHBOR_RADIUS = 10.0      # Radius (Å) to count neighboring Cβ atoms
EXPOSURE_THRESHOLD = 0.5             # Fraction for exposure determination
MAX_NEIGHBORS_BURIED = 24            # Max neighbors for fully buried residue

# Score thresholds for stratification
SCORE_VERY_HIGH_THRESHOLD = 600
SCORE_HIGH_THRESHOLD = 400
SCORE_MEDIUM_THRESHOLD = 200

# Residue classifications
HYDROPHOBIC_RESIDUES = {'A', 'V', 'L', 'I', 'M', 'F', 'C', 'P', 'Y'}
ACIDIC_RESIDUES = {'D', 'E'}

# Parallel processing settings
MAX_WORKERS = 10

# =============================================================================

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Global session for connection pooling
_session = None
_session_lock = Lock()


def get_session() -> requests.Session:
    global _session
    if _session is None:
        with _session_lock:
            if _session is None:
                _session = requests.Session()
                adapter = requests.adapters.HTTPAdapter(
                    pool_connections=MAX_WORKERS,
                    pool_maxsize=MAX_WORKERS * 2
                )
                _session.mount('https://', adapter)
                _session.mount('http://', adapter)
    return _session


AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O'
}


def safe_get(url: str, timeout: int = 30) -> requests.Response | None:
    try:
        r = get_session().get(url, timeout=timeout)
        return r if r.status_code == 200 else None
    except Exception:
        return None


def calc_distance(c1: tuple, c2: tuple) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def calculate_exposure(coords_cb: dict, pos: int, all_cb_coords: list) -> bool:
    """Determine if a residue is surface-exposed based on neighbor counting."""
    if pos not in coords_cb:
        return False

    target_coord = coords_cb[pos]
    neighbor_count = 0

    for cb_coord in all_cb_coords:
        if cb_coord == target_coord:
            continue
        dist = calc_distance(target_coord, cb_coord)
        if dist <= EXPOSURE_NEIGHBOR_RADIUS:
            neighbor_count += 1

    exposure_cutoff = MAX_NEIGHBORS_BURIED * (1 - EXPOSURE_THRESHOLD)
    return neighbor_count < exposure_cutoff


# =============================================================================
# ALPHAFOLD STRUCTURE FETCHING FOR CONTROL ANALYSIS
# =============================================================================

def fetch_alphafold(uid: str) -> dict | None:
    """Fetch AlphaFold structure for a UniProt accession."""
    r = safe_get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}")
    if not r:
        return None
    try:
        data = r.json()
        if not data:
            return None
        cif_url = data[0].get('cifUrl') if data else None
        if not cif_url:
            return None
        cif_r = safe_get(cif_url)
        return parse_cif(cif_r.text) if cif_r else None
    except Exception:
        return None


def parse_cif(content: str) -> dict:
    """Parse mmCIF file to extract coordinates, pLDDT values, and exposure data."""
    plddt, seq, coords_ca, coords_cb, coords_nz = {}, {}, {}, {}, {}
    acidic_coords_sidechain = {}
    acidic_coords_ca = {}
    in_atom, headers, idx = False, [], {}

    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('_atom_site.'):
            in_atom = True
            name = line.split('.')[1].split()[0]
            idx[name] = len(headers)
            headers.append(name)
        elif in_atom and line.startswith(('ATOM', 'HETATM')):
            p = line.split()
            try:
                atom_name = p[idx['label_atom_id']]
                rn = int(p[idx['label_seq_id']])
                aa3 = p[idx['label_comp_id']]
                aa = AA_3TO1.get(aa3, 'X')
                coord = (float(p[idx['Cartn_x']]), float(p[idx['Cartn_y']]), float(p[idx['Cartn_z']]))

                if atom_name == 'CA':
                    plddt[rn] = float(p[idx['B_iso_or_equiv']])
                    seq[rn] = aa
                    coords_ca[rn] = coord
                    if aa3 == 'GLY':
                        coords_cb[rn] = coord
                    if aa3 in ('ASP', 'GLU'):
                        acidic_coords_ca[rn] = coord
                elif atom_name == 'CB':
                    coords_cb[rn] = coord
                elif atom_name == 'NZ' and aa3 == 'LYS':
                    coords_nz[rn] = coord
                elif atom_name == 'CG' and aa3 == 'ASP':
                    acidic_coords_sidechain[rn] = coord
                elif atom_name == 'CD' and aa3 == 'GLU':
                    acidic_coords_sidechain[rn] = coord
            except Exception:
                continue
        elif in_atom and line.startswith('#'):
            break

    lysine_positions = [pos for pos, aa in seq.items() if aa == 'K']
    acidic_positions = sorted(set(acidic_coords_sidechain.keys()) | set(acidic_coords_ca.keys()))

    # Calculate exposure for all residues
    all_cb_coords = list(coords_cb.values())
    exposure = {}
    for pos in coords_cb:
        exposure[pos] = calculate_exposure(coords_cb, pos, all_cb_coords)

    return {
        'plddt': plddt,
        'sequence': seq,
        'coords_ca': coords_ca,
        'coords_cb': coords_cb,
        'coords_nz': coords_nz,
        'acidic_coords_sidechain': acidic_coords_sidechain,
        'acidic_coords_ca': acidic_coords_ca,
        'acidic_positions': acidic_positions,
        'lysine_positions': lysine_positions,
        'exposure': exposure
    }


def analyze_lysine(struct: dict, pos: int) -> dict:
    """Analyze a lysine residue focusing on exposed acidic residues.

    Distance calculation depends on flexibility of BOTH residues:
    - If EITHER lysine OR acidic has pLDDT < threshold: use Cα-Cα distance
    - If BOTH lysine AND acidic have pLDDT >= threshold: use NZ-CG/CD distance
    """
    res = {
        'position': pos,
        'plddt_site': None,
        'plddt_window_avg': None,
        'lys_exposed': None,
        'n_exposed_acidic_within': 0,
        'close_exposed_acidic': None,
        'flexible': None, 'structured': None,
        'Category': None
    }

    if not struct or pos not in struct['plddt']:
        return res

    plddt = struct['plddt']
    seq = struct['sequence']
    coords_ca = struct['coords_ca']
    coords_nz = struct['coords_nz']
    acidic_coords_sidechain = struct['acidic_coords_sidechain']
    acidic_coords_ca = struct['acidic_coords_ca']
    acidic_positions = struct['acidic_positions']
    exposure = struct['exposure']

    res['plddt_site'] = plddt.get(pos)
    lys_plddt = plddt.get(pos, 0)

    half_window = PLDDT_WINDOW_LENGTH // 2
    window_vals = [plddt[i] for i in range(pos - half_window, pos + half_window + 1) if i in plddt]
    res['plddt_window_avg'] = sum(window_vals) / len(window_vals) if window_vals else None

    # Is lysine exposed?
    res['lys_exposed'] = 'Yes' if exposure.get(pos, False) else None

    # Determine flexibility based on window average (for category classification)
    if res['plddt_window_avg'] is not None:
        if res['plddt_window_avg'] < PLDDT_THRESHOLD:
            res['flexible'] = 'Yes'
        else:
            res['structured'] = 'Yes'

    # Calculate distances to EXPOSED acidic residues
    # Method depends on BOTH lysine AND acidic pLDDT
    lys_ca = coords_ca.get(pos)
    lys_nz = coords_nz.get(pos)

    if acidic_positions and (lys_ca or lys_nz):
        n_exposed_within = 0

        for ac_pos in acidic_positions:
            # Only consider EXPOSED acidic residues
            if not exposure.get(ac_pos, False):
                continue

            ac_plddt = plddt.get(ac_pos, 0)

            # Use Cα-Cα if EITHER residue is flexible (pLDDT < threshold)
            # Use NZ-CG/CD only if BOTH are structured (pLDDT >= threshold)
            if lys_plddt < PLDDT_THRESHOLD or ac_plddt < PLDDT_THRESHOLD:
                # At least one is flexible - use Cα-Cα
                if lys_ca is None or ac_pos not in acidic_coords_ca:
                    continue
                d = calc_distance(lys_ca, acidic_coords_ca[ac_pos])
                distance_threshold = DISTANCE_THRESHOLD_FLEXIBLE
            else:
                # Both are structured - use NZ-CG/CD
                if lys_nz is None or ac_pos not in acidic_coords_sidechain:
                    continue
                d = calc_distance(lys_nz, acidic_coords_sidechain[ac_pos])
                distance_threshold = DISTANCE_THRESHOLD_STRUCTURED

            if d <= distance_threshold:
                n_exposed_within += 1

        res['n_exposed_acidic_within'] = n_exposed_within
        if n_exposed_within > 0:
            res['close_exposed_acidic'] = 'Yes'

    # Category assignment
    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"
    suffix = "_exposed" if res['close_exposed_acidic'] == 'Yes' else "_noexposed"
    res['Category'] = prefix + suffix

    return res


def get_all_lysines_analysis(uniprot_ids: list) -> pd.DataFrame:
    """Fetch structures and analyze all lysines for control comparison."""
    unique_ids = list(set(uniprot_ids))
    logger.info(f"Fetching {len(unique_ids)} unique structures for control analysis...")

    all_lysines = []
    completed = 0

    def process_uid(uid):
        struct = fetch_alphafold(uid)
        if not struct:
            return []
        lysines = []
        for pos in struct['lysine_positions']:
            analysis = analyze_lysine(struct, pos)
            analysis['uniprot_id'] = uid
            lysines.append(analysis)
        return lysines

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_uid, uid): uid for uid in unique_ids}
        for future in as_completed(futures):
            lysines = future.result()
            all_lysines.extend(lysines)
            completed += 1
            if completed % 50 == 0 or completed == len(unique_ids):
                logger.info(f"Control analysis: {completed}/{len(unique_ids)} structures processed")

    return pd.DataFrame(all_lysines)


# =============================================================================
# STATISTICAL ANALYSIS FUNCTIONS
# =============================================================================

def count_categories(df: pd.DataFrame) -> dict:
    """Count sites in each category (exposed-focused)."""
    categories = [
        'Flexible_exposed', 'Flexible_noexposed',
        'Structured_exposed', 'Structured_noexposed'
    ]

    cat_col = 'Category' if 'Category' in df.columns else 'category'

    counts = {}
    for cat in categories:
        counts[cat] = len(df[df[cat_col] == cat])

    counts['Total_Flexible'] = counts['Flexible_exposed'] + counts['Flexible_noexposed']
    counts['Total_Structured'] = counts['Structured_exposed'] + counts['Structured_noexposed']
    counts['Total_Exposed'] = counts['Flexible_exposed'] + counts['Structured_exposed']
    counts['Total_Noexposed'] = counts['Flexible_noexposed'] + counts['Structured_noexposed']
    counts['Total'] = counts['Total_Flexible'] + counts['Total_Structured']

    return counts


def calculate_rates(counts: dict) -> dict:
    """Calculate rates within Flexible and Structured."""
    rates = {}
    for prefix in ['Flexible', 'Structured']:
        total = counts[f'Total_{prefix}']
        if total > 0:
            for suffix in ['exposed', 'noexposed']:
                key = f'{prefix}_{suffix}'
                rates[f'{key}_rate'] = counts[key] / total
        else:
            for suffix in ['exposed', 'noexposed']:
                key = f'{prefix}_{suffix}'
                rates[f'{key}_rate'] = None
    return rates


def compare_flex_vs_struct(df: pd.DataFrame, label: str = "") -> list:
    """Chi-square test comparing category distribution between Flexible and Structured."""
    results = []
    results.append({'Metric': f'=== {label} FLEXIBLE vs STRUCTURED ===', 'Value': ''})

    cat_col = 'Category' if 'Category' in df.columns else 'category'
    if cat_col not in df.columns:
        results.append({'Metric': 'Error', 'Value': 'No Category column found'})
        return results

    counts = count_categories(df)
    rates = calculate_rates(counts)

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'Category', 'Value': 'Flexible | Structured'})
    for suffix in ['exposed', 'noexposed']:
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
    flex_counts = [counts['Flexible_exposed'], counts['Flexible_noexposed']]
    struct_counts = [counts['Structured_exposed'], counts['Structured_noexposed']]

    if sum(flex_counts) > 0 and sum(struct_counts) > 0:
        results.append({'Metric': '', 'Value': ''})
        try:
            chi2, p, dof, expected = stats.chi2_contingency([flex_counts, struct_counts])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Chi-square test', 'Value': f'χ²={chi2:.2f}, df={dof}, p={p:.2e} {sig}'})
        except Exception as e:
            results.append({'Metric': 'Chi-square error', 'Value': str(e)})

        # Fisher exact test
        results.append({'Metric': '', 'Value': ''})
        try:
            odds, p = stats.fisher_exact([flex_counts, struct_counts])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Fisher exact (exposed vs noexposed)', 'Value': f'OR={odds:.2f}, p={p:.2e} {sig}'})
        except:
            pass

    return results


def compare_sites_vs_control(sites_df: pd.DataFrame, control_df: pd.DataFrame) -> list:
    """Compare SUMO sites vs all lysines (control) focusing on exposed acidics."""
    results = []
    results.append({'Metric': '=== SUMO SITES vs ALL LYSINES (CONTROL) ===', 'Value': ''})

    sites_counts = count_categories(sites_df)
    control_counts = count_categories(control_df)

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'N sites', 'Value': sites_counts['Total']})
    results.append({'Metric': 'N control lysines', 'Value': control_counts['Total']})

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': '--- Overall category comparison ---', 'Value': ''})
    results.append({'Metric': 'Category', 'Value': 'Sites | Control'})

    for cat in ['Flexible_exposed', 'Flexible_noexposed', 'Structured_exposed', 'Structured_noexposed']:
        s_count = sites_counts[cat]
        c_count = control_counts[cat]
        s_pct = 100 * s_count / sites_counts['Total'] if sites_counts['Total'] > 0 else 0
        c_pct = 100 * c_count / control_counts['Total'] if control_counts['Total'] > 0 else 0
        results.append({'Metric': cat, 'Value': f'{s_count} ({s_pct:.1f}%) | {c_count} ({c_pct:.1f}%)'})

    # Chi-square comparing overall distributions
    results.append({'Metric': '', 'Value': ''})
    sites_all = [sites_counts[c] for c in ['Flexible_exposed', 'Flexible_noexposed', 'Structured_exposed', 'Structured_noexposed']]
    control_all = [control_counts[c] for c in ['Flexible_exposed', 'Flexible_noexposed', 'Structured_exposed', 'Structured_noexposed']]

    if sum(sites_all) > 0 and sum(control_all) > 0:
        try:
            chi2, p, dof, expected = stats.chi2_contingency([sites_all, control_all])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Chi-square (all categories)', 'Value': f'χ²={chi2:.2f}, df={dof}, p={p:.2e} {sig}'})
        except Exception as e:
            results.append({'Metric': 'Chi-square error', 'Value': str(e)})

    # Compare specific features
    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': '--- Feature enrichment in SUMO sites vs control ---', 'Value': ''})

    features = [
        ('Flexible', 'flexible', 'Flexible'),
        ('Lys_exposed', 'lys_exposed', 'Lys_exposed'),
        ('Close_exposed_acidic', 'close_exposed_acidic', 'Close_exposed_acidic'),
    ]

    for name, ctrl_col, site_col in features:
        try:
            s_col = site_col if site_col in sites_df.columns else ctrl_col
            c_col = ctrl_col if ctrl_col in control_df.columns else site_col

            s_yes = (sites_df[s_col] == 'Yes').sum() if s_col in sites_df.columns else 0
            c_yes = (control_df[c_col] == 'Yes').sum() if c_col in control_df.columns else 0

            s_no = len(sites_df) - s_yes
            c_no = len(control_df) - c_yes

            s_rate = s_yes / len(sites_df) if len(sites_df) > 0 else 0
            c_rate = c_yes / len(control_df) if len(control_df) > 0 else 0

            odds, p = stats.fisher_exact([[s_yes, s_no], [c_yes, c_no]])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({
                'Metric': name,
                'Value': f'Sites={s_rate:.1%}, Control={c_rate:.1%}, OR={odds:.2f}, p={p:.2e} {sig}'
            })
        except Exception as e:
            results.append({'Metric': name, 'Value': f'Error: {e}'})

    return results


def analyze_score_predictors(df: pd.DataFrame) -> list:
    """Analyze which features predict high score (exposed acidic focus)."""
    results = []
    results.append({'Metric': '=== PREDICTORS OF HIGH SCORE (Exposed Acidic Focus) ===', 'Value': ''})
    results.append({'Metric': f'High score defined as >= {SCORE_HIGH_THRESHOLD}', 'Value': ''})
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
    results.append({'Metric': 'High score sites', 'Value': f'{n_high} ({100*n_high/len(valid):.1f}%)'})
    results.append({'Metric': 'Low score sites', 'Value': f'{n_low} ({100*n_low/len(valid):.1f}%)'})
    results.append({'Metric': '', 'Value': ''})

    # Define predictors (exposed acidic focus)
    predictors = []
    if 'Flexible' in valid.columns:
        predictors.append(('Flexible', valid['Flexible'] == 'Yes'))
    if 'Structured' in valid.columns:
        predictors.append(('Structured', valid['Structured'] == 'Yes'))
    if 'Lys_exposed' in valid.columns:
        predictors.append(('Lys_exposed', valid['Lys_exposed'] == 'Yes'))
    if 'Close_exposed_acidic' in valid.columns:
        predictors.append(('Close_exposed_acidic', valid['Close_exposed_acidic'] == 'Yes'))

    results.append({'Metric': '--- Binary Predictors (Odds Ratios) ---', 'Value': ''})
    results.append({'Metric': 'Predictor', 'Value': 'Rate_High | Rate_Low | OR | p-value'})

    for name, mask in predictors:
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
                results.append({
                    'Metric': name,
                    'Value': f'{rate_high:.1%} | {rate_low:.1%} | {odds:.2f} | N/A'
                })
        else:
            results.append({
                'Metric': name,
                'Value': f'{rate_high:.1%} | {rate_low:.1%} | N/A | N/A'
            })

    # N exposed acidics as predictor
    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': '--- N_exposed_acidic_within as predictor ---', 'Value': ''})

    n_acidic_col = 'N_exposed_acidic_within'
    if n_acidic_col in valid.columns:
        valid['n_exposed_acidic'] = pd.to_numeric(valid[n_acidic_col], errors='coerce').fillna(0).astype(int)

        # Correlation with score
        rho, p = stats.spearmanr(valid['n_exposed_acidic'], valid[score_col])
        results.append({'Metric': 'Spearman: N_exposed_acidic vs Score', 'Value': f'ρ={rho:.3f}, p={p:.2e}'})

        # Point-biserial with high_score
        r, p = stats.pointbiserialr(valid['high_score'], valid['n_exposed_acidic'])
        results.append({'Metric': 'Point-biserial: N_exposed_acidic vs high_score', 'Value': f'r={r:.3f}, p={p:.2e}'})

        # Mean N_exposed_acidic by score group
        results.append({'Metric': '', 'Value': ''})
        high_mean = valid[valid['high_score'] == 1]['n_exposed_acidic'].mean()
        low_mean = valid[valid['high_score'] == 0]['n_exposed_acidic'].mean()
        t_stat, t_p = stats.ttest_ind(
            valid[valid['high_score'] == 1]['n_exposed_acidic'],
            valid[valid['high_score'] == 0]['n_exposed_acidic']
        )
        results.append({'Metric': 'Mean N_exposed_acidic (high score)', 'Value': f'{high_mean:.2f}'})
        results.append({'Metric': 'Mean N_exposed_acidic (low score)', 'Value': f'{low_mean:.2f}'})
        results.append({'Metric': 't-test', 'Value': f't={t_stat:.2f}, p={t_p:.2e}'})

        # Test 2+ exposed acidics
        results.append({'Metric': '', 'Value': ''})
        has_2plus = (valid['n_exposed_acidic'] >= 2).astype(int)
        tp = ((has_2plus == 1) & (valid['high_score'] == 1)).sum()
        fp = ((has_2plus == 1) & (valid['high_score'] == 0)).sum()
        fn = ((has_2plus == 0) & (valid['high_score'] == 1)).sum()
        tn = ((has_2plus == 0) & (valid['high_score'] == 0)).sum()

        if fp > 0 and fn > 0:
            odds = (tp * tn) / (fp * fn)
            _, p = stats.fisher_exact([[tp, fp], [fn, tn]])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Has 2+ exposed acidics within threshold', 'Value': f'OR={odds:.2f}, p={p:.2e} {sig}'})

    # ROC-style metrics
    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': '--- Classifier Performance Metrics ---', 'Value': ''})
    results.append({'Metric': 'Predictor', 'Value': 'Sens | Spec | PPV | NPV | Accuracy'})

    for name, mask in predictors:
        pred = mask.astype(int)
        tp = ((pred == 1) & (valid['high_score'] == 1)).sum()
        fp = ((pred == 1) & (valid['high_score'] == 0)).sum()
        fn = ((pred == 0) & (valid['high_score'] == 1)).sum()
        tn = ((pred == 0) & (valid['high_score'] == 0)).sum()

        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        acc = (tp + tn) / len(valid)

        results.append({
            'Metric': name,
            'Value': f'{sens:.2f} | {spec:.2f} | {ppv:.2f} | {npv:.2f} | {acc:.2f}'
        })

    return results


def analyze_score_predictors_by_flexibility(df: pd.DataFrame) -> dict:
    """Analyze which features predict high score, separately for Flexible and Structured sites.

    Returns a dict with two DataFrames: one for Flexible sites, one for Structured sites.
    """
    sheets = {}

    score_col = 'Score (SUMO site)'
    if score_col not in df.columns:
        return sheets

    for flex_type, flex_filter in [('Flexible', df['Flexible'] == 'Yes'),
                                    ('Structured', df['Structured'] == 'Yes')]:
        results = []
        results.append({'Metric': f'=== PREDICTORS OF HIGH SCORE ({flex_type.upper()} SITES ONLY) ===', 'Value': ''})
        results.append({'Metric': f'High score defined as >= {SCORE_HIGH_THRESHOLD}', 'Value': ''})
        results.append({'Metric': '', 'Value': ''})

        subset = df[flex_filter & df['pLDDT_site'].notna() & df[score_col].notna()].copy()
        subset[score_col] = pd.to_numeric(subset[score_col], errors='coerce')
        subset = subset[subset[score_col].notna()]

        if len(subset) < 20:
            results.append({'Metric': 'Error', 'Value': f'Too few {flex_type} samples ({len(subset)})'})
            sheets[f'Predictors_{flex_type}'] = pd.DataFrame(results)
            continue

        subset['high_score'] = (subset[score_col] >= SCORE_HIGH_THRESHOLD).astype(int)
        n_high = subset['high_score'].sum()
        n_low = len(subset) - n_high

        results.append({'Metric': f'Total {flex_type} sites', 'Value': len(subset)})
        results.append({'Metric': 'High score sites', 'Value': f'{n_high} ({100*n_high/len(subset):.1f}%)'})
        results.append({'Metric': 'Low score sites', 'Value': f'{n_low} ({100*n_low/len(subset):.1f}%)'})
        results.append({'Metric': '', 'Value': ''})

        # Define predictors (exposed acidic focus, excluding Flexible/Structured)
        predictors = []
        if 'Lys_exposed' in subset.columns:
            predictors.append(('Lys_exposed', subset['Lys_exposed'] == 'Yes'))
        if 'Close_exposed_acidic' in subset.columns:
            predictors.append(('Close_exposed_acidic', subset['Close_exposed_acidic'] == 'Yes'))

        results.append({'Metric': '--- Binary Predictors (Odds Ratios) ---', 'Value': ''})
        results.append({'Metric': 'Predictor', 'Value': 'Rate_High | Rate_Low | OR | p-value'})

        for name, mask in predictors:
            pred = mask.astype(int)
            tp = ((pred == 1) & (subset['high_score'] == 1)).sum()
            fp = ((pred == 1) & (subset['high_score'] == 0)).sum()
            fn = ((pred == 0) & (subset['high_score'] == 1)).sum()
            tn = ((pred == 0) & (subset['high_score'] == 0)).sum()

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
                    results.append({
                        'Metric': name,
                        'Value': f'{rate_high:.1%} | {rate_low:.1%} | {odds:.2f} | N/A'
                    })
            else:
                results.append({
                    'Metric': name,
                    'Value': f'{rate_high:.1%} | {rate_low:.1%} | N/A | N/A'
                })

        # N exposed acidic analysis
        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': '--- N_exposed_acidic_within as predictor ---', 'Value': ''})

        n_acidic_col = 'N_exposed_acidic_within'
        if n_acidic_col in subset.columns:
            subset['n_exposed_acidic'] = pd.to_numeric(subset[n_acidic_col], errors='coerce').fillna(0).astype(int)

            # Correlation with score
            if len(subset) >= 10:
                rho, p = stats.spearmanr(subset['n_exposed_acidic'], subset[score_col])
                results.append({'Metric': 'Spearman: N_exposed_acidic vs Score', 'Value': f'ρ={rho:.3f}, p={p:.2e}'})

                # Point-biserial with high_score
                r, p = stats.pointbiserialr(subset['high_score'], subset['n_exposed_acidic'])
                results.append({'Metric': 'Point-biserial: N_exposed_acidic vs high_score', 'Value': f'r={r:.3f}, p={p:.2e}'})

                # Mean N_exposed_acidic by score group
                results.append({'Metric': '', 'Value': ''})
                high_mean = subset[subset['high_score'] == 1]['n_exposed_acidic'].mean()
                low_mean = subset[subset['high_score'] == 0]['n_exposed_acidic'].mean()

                if n_high >= 2 and n_low >= 2:
                    t_stat, t_p = stats.ttest_ind(
                        subset[subset['high_score'] == 1]['n_exposed_acidic'],
                        subset[subset['high_score'] == 0]['n_exposed_acidic']
                    )
                    results.append({'Metric': 'Mean N_exposed_acidic (high score)', 'Value': f'{high_mean:.2f}'})
                    results.append({'Metric': 'Mean N_exposed_acidic (low score)', 'Value': f'{low_mean:.2f}'})
                    results.append({'Metric': 't-test', 'Value': f't={t_stat:.2f}, p={t_p:.2e}'})

        # Classifier performance metrics
        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': '--- Classifier Performance Metrics ---', 'Value': ''})
        results.append({'Metric': 'Predictor', 'Value': 'Sens | Spec | PPV | NPV | Accuracy'})

        for name, mask in predictors:
            pred = mask.astype(int)
            tp = ((pred == 1) & (subset['high_score'] == 1)).sum()
            fp = ((pred == 1) & (subset['high_score'] == 0)).sum()
            fn = ((pred == 0) & (subset['high_score'] == 1)).sum()
            tn = ((pred == 0) & (subset['high_score'] == 0)).sum()

            sens = tp / (tp + fn) if (tp + fn) > 0 else 0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0
            ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
            npv = tn / (tn + fn) if (tn + fn) > 0 else 0
            acc = (tp + tn) / len(subset) if len(subset) > 0 else 0

            results.append({
                'Metric': name,
                'Value': f'{sens:.2f} | {spec:.2f} | {ppv:.2f} | {npv:.2f} | {acc:.2f}'
            })

        sheets[f'Predictors_{flex_type}'] = pd.DataFrame(results)

    return sheets


def perform_global_analysis(df: pd.DataFrame, control_df: pd.DataFrame = None) -> pd.DataFrame:
    """Perform global analysis on all sites (exposed acidic focus)."""
    results = []

    results.append({'Metric': '=== GLOBAL ANALYSIS (Exposed Acidic Focus) ===', 'Value': ''})
    results.append({'Metric': '', 'Value': ''})

    counts = count_categories(df)
    rates = calculate_rates(counts)

    results.append({'Metric': '--- SUMO Sites Category Counts ---', 'Value': ''})
    for cat in ['Flexible_exposed', 'Flexible_noexposed']:
        rate = rates.get(f'{cat}_rate', 0) or 0
        results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%} of Flexible)'})

    results.append({'Metric': '', 'Value': ''})
    for cat in ['Structured_exposed', 'Structured_noexposed']:
        rate = rates.get(f'{cat}_rate', 0) or 0
        results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%} of Structured)'})

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'Total Flexible', 'Value': counts['Total_Flexible']})
    results.append({'Metric': 'Total Structured', 'Value': counts['Total_Structured']})
    results.append({'Metric': 'Total with exposed acidic', 'Value': counts['Total_Exposed']})
    results.append({'Metric': 'Total', 'Value': counts['Total']})

    # Flex vs Struct comparison
    results.append({'Metric': '', 'Value': ''})
    results.extend(compare_flex_vs_struct(df, "SITES:"))

    # Control analysis
    if control_df is not None and len(control_df) > 0:
        cat_col = 'Category' if 'Category' in control_df.columns else 'category'
        if cat_col not in control_df.columns:
            results.append({'Metric': '', 'Value': ''})
            results.append({'Metric': 'Control Error', 'Value': 'No category data in control'})
            return pd.DataFrame(results)

        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': '', 'Value': ''})

        control_counts = count_categories(control_df)
        control_rates = calculate_rates(control_counts)

        results.append({'Metric': '--- Control (All Lysines) Category Counts ---', 'Value': ''})
        for cat in ['Flexible_exposed', 'Flexible_noexposed']:
            rate = control_rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{control_counts[cat]} ({rate:.1%} of Flexible)'})

        results.append({'Metric': '', 'Value': ''})
        for cat in ['Structured_exposed', 'Structured_noexposed']:
            rate = control_rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{control_counts[cat]} ({rate:.1%} of Structured)'})

        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': 'Total Flexible', 'Value': control_counts['Total_Flexible']})
        results.append({'Metric': 'Total Structured', 'Value': control_counts['Total_Structured']})
        results.append({'Metric': 'Total with exposed acidic', 'Value': control_counts['Total_Exposed']})
        results.append({'Metric': 'Total', 'Value': control_counts['Total']})

        # Flex vs Struct for control
        results.append({'Metric': '', 'Value': ''})
        results.extend(compare_flex_vs_struct(control_df, "CONTROL:"))

        # Sites vs Control comparison
        results.append({'Metric': '', 'Value': ''})
        results.extend(compare_sites_vs_control(df, control_df))

    return pd.DataFrame(results)


def perform_stratified_analysis(df: pd.DataFrame, score_col: str = 'Score (SUMO site)') -> dict:
    """Perform analysis stratified by score groups (exposed acidic focus)."""
    sheets = {}

    if score_col not in df.columns:
        return sheets

    df[score_col] = pd.to_numeric(df[score_col], errors='coerce')

    strata = {
        'VeryHigh': df[df[score_col] >= SCORE_VERY_HIGH_THRESHOLD],
        'High': df[(df[score_col] >= SCORE_HIGH_THRESHOLD) & (df[score_col] < SCORE_VERY_HIGH_THRESHOLD)],
        'Medium': df[(df[score_col] >= SCORE_MEDIUM_THRESHOLD) & (df[score_col] < SCORE_HIGH_THRESHOLD)],
        'Low': df[df[score_col] < SCORE_MEDIUM_THRESHOLD]
    }

    for name, subset in strata.items():
        if len(subset) == 0:
            continue

        results = []
        results.append({'Metric': f'=== {name.upper()} SCORE ANALYSIS ===', 'Value': ''})
        results.append({'Metric': 'N sites', 'Value': len(subset)})
        results.append({'Metric': '', 'Value': ''})

        counts = count_categories(subset)
        rates = calculate_rates(counts)

        results.append({'Metric': '--- Category Counts ---', 'Value': ''})
        for cat in ['Flexible_exposed', 'Flexible_noexposed']:
            rate = rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%})'})

        results.append({'Metric': '', 'Value': ''})
        for cat in ['Structured_exposed', 'Structured_noexposed']:
            rate = rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%})'})

        results.append({'Metric': '', 'Value': ''})
        results.extend(compare_flex_vs_struct(subset, name.upper()))

        sheets[f'Score_{name}'] = pd.DataFrame(results)

    return sheets


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def analyze_file(input_file: str, output_file: str = None):
    """Load core output and perform all analyses (exposed acidic focus)."""
    print(f"\n{'=' * 60}\nSUMO Site Analyzer - Statistical Analysis (Exposed Acidic Focus)\n{'=' * 60}\n")
    print(f"Parameters:")
    print(f"  pLDDT threshold: {PLDDT_THRESHOLD}")
    print(f"  Distance threshold (flexible, Cα-Cα): {DISTANCE_THRESHOLD_FLEXIBLE} Å")
    print(f"  Distance threshold (structured, NZ-CG/CD): {DISTANCE_THRESHOLD_STRUCTURED} Å")
    print(f"  Exposure neighbor radius: {EXPOSURE_NEIGHBOR_RADIUS} Å")
    print(f"  Exposure threshold: {EXPOSURE_THRESHOLD} ({int(EXPOSURE_THRESHOLD*100)}%)")
    print(f"  Score thresholds: {SCORE_VERY_HIGH_THRESHOLD}/{SCORE_HIGH_THRESHOLD}/{SCORE_MEDIUM_THRESHOLD}")
    print()

    # Load data
    xls = pd.ExcelFile(input_file)
    sheet_name = 'Sites' if 'Sites' in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(input_file, sheet_name=sheet_name)

    print(f"Loaded {len(df)} sites from {input_file}")

    # Get valid sites (those with pLDDT data)
    valid_df = df[df['pLDDT_site'].notna()].copy()
    print(f"Valid sites with structure data: {len(valid_df)}")

    # Get unique UniProt IDs for control analysis
    uid_col = 'UniProt_ID_used'
    if uid_col in valid_df.columns:
        uniprot_ids = valid_df[uid_col].dropna().unique().tolist()
        print(f"Unique proteins: {len(uniprot_ids)}")

        # Fetch control data (all lysines)
        print("\nFetching control data (all lysines in analyzed structures)...")
        control_df = get_all_lysines_analysis(uniprot_ids)
        print(f"Control lysines analyzed: {len(control_df)}")
    else:
        control_df = pd.DataFrame()
        print("Warning: UniProt_ID_used column not found, skipping control analysis")

    # Prepare output
    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_analysis_exposed.xlsx'

    print("\nGenerating analysis sheets...")

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Global analysis
        global_df = perform_global_analysis(valid_df, control_df if len(control_df) > 0 else None)
        global_df.to_excel(writer, sheet_name='Global_Analysis', index=False)

        # Score-stratified analysis
        strata_sheets = perform_stratified_analysis(valid_df)
        for sheet_name, sheet_df in strata_sheets.items():
            sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)

        # Predictors of high score
        predictors_results = analyze_score_predictors(valid_df)
        pd.DataFrame(predictors_results).to_excel(writer, sheet_name='Score_Predictors', index=False)

        # Predictors of high score by flexibility (Flexible and Structured separately)
        flex_struct_predictors = analyze_score_predictors_by_flexibility(valid_df)
        for sheet_name, sheet_df in flex_struct_predictors.items():
            sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)

        # Control comparison sheet (detailed)
        if len(control_df) > 0:
            control_comparison = compare_sites_vs_control(valid_df, control_df)
            pd.DataFrame(control_comparison).to_excel(writer, sheet_name='Sites_vs_Control', index=False)

            # Control Flex vs Struct
            control_flex_struct = compare_flex_vs_struct(control_df, "CONTROL LYSINES:")
            pd.DataFrame(control_flex_struct).to_excel(writer, sheet_name='Control_FlexVsStruct', index=False)

    print(f"\nAnalysis complete. Output: {output_file}")
    return output_file


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nConfigurable parameters (must match sumo_site_core_exposed.py):")
        print(f"  PLDDT_THRESHOLD = {PLDDT_THRESHOLD}")
        print(f"  DISTANCE_THRESHOLD_FLEXIBLE = {DISTANCE_THRESHOLD_FLEXIBLE} (Cα-Cα)")
        print(f"  DISTANCE_THRESHOLD_STRUCTURED = {DISTANCE_THRESHOLD_STRUCTURED} (NZ-CG/CD)")
        print(f"  EXPOSURE_NEIGHBOR_RADIUS = {EXPOSURE_NEIGHBOR_RADIUS}")
        print(f"  EXPOSURE_THRESHOLD = {EXPOSURE_THRESHOLD}")
        print(f"  SCORE_VERY_HIGH_THRESHOLD = {SCORE_VERY_HIGH_THRESHOLD}")
        print(f"  SCORE_HIGH_THRESHOLD = {SCORE_HIGH_THRESHOLD}")
        print(f"  SCORE_MEDIUM_THRESHOLD = {SCORE_MEDIUM_THRESHOLD}")
        sys.exit(1)
    analyze_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
