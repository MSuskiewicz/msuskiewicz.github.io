#!/usr/bin/env python3
"""
SUMO Site Analyzer - Statistical Analysis
==========================================
Analyzes output from sumo_site_core.py with statistical comparisons.

Analyses:
1. Global category counts (Flexible/Structured × consensus/flankacidic/space/none)
2. Score-stratified analysis (VeryHigh, High, Medium, Low)
3. Control analysis: all lysines in analyzed AlphaFold models
4. Statistical comparisons: Flexible vs Structured, Sites vs Control
5. Predictors of high score

Usage:
    python sumo_site_analysis.py <core_output_excel> [output_file]
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
# CONFIGURABLE PARAMETERS (must match sumo_site_core.py)
# =============================================================================

PLDDT_THRESHOLD = 65.0          # pLDDT below this = Flexible, above = Structured
PLDDT_WINDOW_LENGTH = 11        # Window size for pLDDT averaging (centered on site)

# Distance thresholds for ACIDIC residues (Cα-Cα based, minimum and maximum)
DISTANCE_THRESHOLD_MIN = 4.9    # Minimum Cα-Cα distance in Angstroms
DISTANCE_THRESHOLD_MAX = 8.0    # Maximum Cα-Cα distance in Angstroms

# Distance thresholds for HYDROPHOBIC residues (Cα-Cα based)
HYDROPHOBIC_DISTANCE_MIN = 3.0  # Minimum Cα-Cα distance in Angstroms
HYDROPHOBIC_DISTANCE_MAX = 4.5  # Maximum Cα-Cα distance in Angstroms

# Surface exposure parameters (neighbor counting method)
EXPOSURE_DISTANCE_THRESHOLD = 10.0  # Angstroms - radius for counting CB neighbors
EXPOSURE_NEIGHBOR_THRESHOLD = 18    # Fewer than this many neighbors = exposed

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
    """Parse mmCIF file to extract coordinates and pLDDT values.

    Extracts:
    - CA coordinates for all residues (for Cα-Cα distance calculations)
    - CB coordinates for all residues (for exposure calculation)
    - pLDDT values from B-factor column
    - Acidic residue (D, E) coordinates
    - Hydrophobic residue coordinates
    - Secondary structure from _struct_conf category
    """
    plddt, seq, coords_ca, coords_cb = {}, {}, {}, {}
    acidic_coords_ca = {}  # CA for acidic residues
    hydrophobic_coords_ca = {}  # CA for hydrophobic residues
    in_atom, headers, idx = False, [], {}

    # Secondary structure storage
    sec_struct = {}  # residue_num -> 'helix' | 'strand' | 'coil'

    lines = content.split('\n')

    # First pass: extract secondary structure from _struct_conf
    in_struct_conf = False
    struct_conf_headers = []
    struct_conf_idx = {}

    for line in lines:
        line_stripped = line.strip()
        if line_stripped.startswith('_struct_conf.'):
            in_struct_conf = True
            name = line_stripped.split('.')[1].split()[0]
            struct_conf_idx[name] = len(struct_conf_headers)
            struct_conf_headers.append(name)
        elif in_struct_conf and line_stripped and not line_stripped.startswith(('_', '#', 'loop_')):
            if 'conf_type_id' in struct_conf_idx:
                try:
                    parts = line_stripped.split()
                    conf_type = parts[struct_conf_idx['conf_type_id']]
                    beg_seq = int(parts[struct_conf_idx['beg_label_seq_id']])
                    end_seq = int(parts[struct_conf_idx['end_label_seq_id']])

                    ss_type = 'coil'
                    if 'HELX' in conf_type or 'helix' in conf_type.lower():
                        ss_type = 'helix'
                    elif 'STRN' in conf_type or 'strand' in conf_type.lower():
                        ss_type = 'strand'

                    for pos in range(beg_seq, end_seq + 1):
                        sec_struct[pos] = ss_type
                except Exception:
                    pass
        elif in_struct_conf and (line_stripped.startswith('#') or line_stripped.startswith('loop_') or line_stripped.startswith('_')):
            if not line_stripped.startswith('_struct_conf.'):
                in_struct_conf = False

    # Second pass: extract atom coordinates
    for line in lines:
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
                    if aa3 in ('ASP', 'GLU'):
                        acidic_coords_ca[rn] = coord
                    if aa in HYDROPHOBIC_RESIDUES:
                        hydrophobic_coords_ca[rn] = coord
                elif atom_name == 'CB':
                    coords_cb[rn] = coord
            except Exception:
                continue
        elif in_atom and line.startswith('#'):
            break

    lysine_positions = [pos for pos, aa in seq.items() if aa == 'K']
    acidic_positions = sorted(acidic_coords_ca.keys())
    hydrophobic_positions = sorted(hydrophobic_coords_ca.keys())

    # For glycines (no CB), use CA as fallback
    for rn in coords_ca:
        if rn not in coords_cb:
            coords_cb[rn] = coords_ca[rn]

    return {
        'plddt': plddt,
        'sequence': seq,
        'coords_ca': coords_ca,
        'coords_cb': coords_cb,
        'acidic_coords_ca': acidic_coords_ca,
        'acidic_positions': acidic_positions,
        'hydrophobic_coords_ca': hydrophobic_coords_ca,
        'hydrophobic_positions': hydrophobic_positions,
        'lysine_positions': lysine_positions,
        'secondary_structure': sec_struct
    }


def calculate_exposure(struct: dict, pos: int) -> bool:
    """Determine if a residue is surface-exposed using neighbor counting."""
    coords_cb = struct.get('coords_cb', {})
    if pos not in coords_cb:
        return False

    target_cb = coords_cb[pos]
    neighbor_count = 0

    for other_pos, other_cb in coords_cb.items():
        if other_pos == pos:
            continue
        d = calc_distance(target_cb, other_cb)
        if d <= EXPOSURE_DISTANCE_THRESHOLD:
            neighbor_count += 1

    return neighbor_count < EXPOSURE_NEIGHBOR_THRESHOLD


def analyze_lysine(struct: dict, pos: int) -> dict:
    """Analyze a lysine residue (same logic as SUMO site analysis).

    Uses Cα-Cα distance for all calculations with min/max thresholds.
    """
    res = {
        'position': pos,
        'plddt_site': None,
        'plddt_window_avg': None,
        'aa_m2': None, 'aa_m1': None, 'aa_p1': None, 'aa_p2': None,
        'flexible': None, 'structured': None,
        'secondary_structure': None,
        'forward_consensus': None, 'inverse_consensus': None,
        'any_acidic_within_threshold': None,
        'exposed_acidic_within_threshold': None,
        'acidic_in_pm2': None,
        'exposed_acidic_in_pm2': None,
        'any_hydrophobic_within_threshold': None,
        'exposed_hydrophobic_within_threshold': None,
        'Category': None  # Use uppercase to match sites
    }

    if not struct or pos not in struct['plddt']:
        return res

    plddt = struct['plddt']
    seq = struct['sequence']
    coords_ca = struct['coords_ca']
    acidic_coords_ca = struct['acidic_coords_ca']
    acidic_positions = struct['acidic_positions']
    hydrophobic_coords_ca = struct.get('hydrophobic_coords_ca', {})
    hydrophobic_positions = struct.get('hydrophobic_positions', [])
    sec_struct = struct.get('secondary_structure', {})

    res['plddt_site'] = plddt.get(pos)

    half_window = PLDDT_WINDOW_LENGTH // 2
    window_vals = [plddt[i] for i in range(pos - half_window, pos + half_window + 1) if i in plddt]
    res['plddt_window_avg'] = sum(window_vals) / len(window_vals) if window_vals else None

    res['aa_m2'] = seq.get(pos - 2)
    res['aa_m1'] = seq.get(pos - 1)
    res['aa_p1'] = seq.get(pos + 1)
    res['aa_p2'] = seq.get(pos + 2)

    # Determine flexibility based on window average
    if res['plddt_window_avg'] is not None:
        if res['plddt_window_avg'] < PLDDT_THRESHOLD:
            res['flexible'] = 'Yes'
        else:
            res['structured'] = 'Yes'

    # Secondary structure
    ss = sec_struct.get(pos)
    res['secondary_structure'] = ss if ss else 'coil'

    # Calculate Cα-Cα distances to all acidic residues
    lys_ca = coords_ca.get(pos)

    if acidic_positions and lys_ca:
        has_any_acidic = False
        has_exposed_acidic = False
        has_in_pm2 = False
        has_exposed_in_pm2 = False

        for ac_pos in acidic_positions:
            if ac_pos not in acidic_coords_ca:
                continue

            d = calc_distance(lys_ca, acidic_coords_ca[ac_pos])

            # Check if within threshold range (min <= d <= max)
            if DISTANCE_THRESHOLD_MIN <= d <= DISTANCE_THRESHOLD_MAX:
                has_any_acidic = True
                is_exposed = calculate_exposure(struct, ac_pos)
                if is_exposed:
                    has_exposed_acidic = True

                # Check relative position for +/-2
                rel_pos = ac_pos - pos
                if -2 <= rel_pos <= 2:
                    has_in_pm2 = True
                    if is_exposed:
                        has_exposed_in_pm2 = True

        if has_any_acidic:
            res['any_acidic_within_threshold'] = 'Yes'
        if has_exposed_acidic:
            res['exposed_acidic_within_threshold'] = 'Yes'
        if has_in_pm2:
            res['acidic_in_pm2'] = 'Yes'
        if has_exposed_in_pm2:
            res['exposed_acidic_in_pm2'] = 'Yes'

    # Calculate Cα-Cα distances to all hydrophobic residues
    if hydrophobic_positions and lys_ca:
        has_any_hp = False
        has_exposed_hp = False

        for hp_pos in hydrophobic_positions:
            if hp_pos not in hydrophobic_coords_ca:
                continue

            d = calc_distance(lys_ca, hydrophobic_coords_ca[hp_pos])

            # Check if within hydrophobic threshold range
            if HYDROPHOBIC_DISTANCE_MIN <= d <= HYDROPHOBIC_DISTANCE_MAX:
                has_any_hp = True
                if calculate_exposure(struct, hp_pos):
                    has_exposed_hp = True

        if has_any_hp:
            res['any_hydrophobic_within_threshold'] = 'Yes'
        if has_exposed_hp:
            res['exposed_hydrophobic_within_threshold'] = 'Yes'

    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    res['inverse_consensus'] = 'Yes' if res['aa_m2'] in ACIDIC_RESIDUES and res['aa_p1'] in HYDROPHOBIC_RESIDUES else None

    is_consensus = (res['forward_consensus'] == 'Yes' or res['inverse_consensus'] == 'Yes')
    has_any_acidic = (res['any_acidic_within_threshold'] == 'Yes')
    has_exposed_acidic = (res['exposed_acidic_within_threshold'] == 'Yes')

    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"

    if is_consensus:
        suffix = "_consensus"
    elif has_exposed_acidic:
        suffix = "_exposed_acidic"
    elif has_any_acidic:
        suffix = "_buried_acidic"
    else:
        suffix = "_no_acidic"

    res['Category'] = prefix + suffix  # Use uppercase to match sites DataFrame
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
    """Count sites in each category."""
    categories = [
        'Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic',
        'Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic'
    ]

    # Handle both uppercase (from core output) and lowercase (from control analysis) column names
    cat_col = 'Category' if 'Category' in df.columns else 'category'

    counts = {}
    for cat in categories:
        counts[cat] = len(df[df[cat_col] == cat])

    counts['Total_Flexible'] = sum(counts[c] for c in categories[:4])
    counts['Total_Structured'] = sum(counts[c] for c in categories[4:])
    counts['Total'] = counts['Total_Flexible'] + counts['Total_Structured']

    return counts


def calculate_rates(counts: dict) -> dict:
    """Calculate rates within Flexible and Structured."""
    rates = {}
    for prefix in ['Flexible', 'Structured']:
        total = counts[f'Total_{prefix}']
        if total > 0:
            for suffix in ['consensus', 'exposed_acidic', 'buried_acidic', 'no_acidic']:
                key = f'{prefix}_{suffix}'
                rates[f'{key}_rate'] = counts[key] / total
        else:
            for suffix in ['consensus', 'exposed_acidic', 'buried_acidic', 'no_acidic']:
                key = f'{prefix}_{suffix}'
                rates[f'{key}_rate'] = None
    return rates


def compare_flex_vs_struct(df: pd.DataFrame, label: str = "") -> list:
    """Chi-square test comparing category distribution between Flexible and Structured."""
    results = []
    results.append({'Metric': f'=== {label} FLEXIBLE vs STRUCTURED ===', 'Value': ''})

    # Handle case where no valid category data exists
    cat_col = 'Category' if 'Category' in df.columns else 'category'
    if cat_col not in df.columns:
        results.append({'Metric': 'Error', 'Value': 'No Category column found'})
        return results

    counts = count_categories(df)
    rates = calculate_rates(counts)

    # Show counts
    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'Category', 'Value': 'Flexible | Structured'})
    for suffix in ['consensus', 'exposed_acidic', 'buried_acidic', 'no_acidic']:
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
    flex_counts = [counts[f'Flexible_{s}'] for s in ['consensus', 'exposed_acidic', 'buried_acidic', 'no_acidic']]
    struct_counts = [counts[f'Structured_{s}'] for s in ['consensus', 'exposed_acidic', 'buried_acidic', 'no_acidic']]

    if sum(flex_counts) > 0 and sum(struct_counts) > 0:
        results.append({'Metric': '', 'Value': ''})
        try:
            chi2, p, dof, expected = stats.chi2_contingency([flex_counts, struct_counts])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Chi-square test', 'Value': f'χ²={chi2:.2f}, df={dof}, p={p:.2e} {sig}'})
        except Exception as e:
            results.append({'Metric': 'Chi-square error', 'Value': str(e)})

        # Fisher exact for each category (Flexible vs Structured)
        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': '--- Fisher exact tests (each category) ---', 'Value': ''})
        for i, suffix in enumerate(['consensus', 'exposed_acidic', 'buried_acidic', 'no_acidic']):
            f_yes, f_no = flex_counts[i], sum(flex_counts) - flex_counts[i]
            s_yes, s_no = struct_counts[i], sum(struct_counts) - struct_counts[i]
            try:
                odds, p = stats.fisher_exact([[f_yes, f_no], [s_yes, s_no]])
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                results.append({'Metric': suffix, 'Value': f'OR={odds:.2f}, p={p:.2e} {sig}'})
            except:
                pass

    return results


def compare_sites_vs_control(sites_df: pd.DataFrame, control_df: pd.DataFrame) -> list:
    """Compare SUMO sites vs all lysines (control)."""
    results = []
    results.append({'Metric': '=== SUMO SITES vs ALL LYSINES (CONTROL) ===', 'Value': ''})

    sites_counts = count_categories(sites_df)
    control_counts = count_categories(control_df)

    sites_rates = calculate_rates(sites_counts)
    control_rates = calculate_rates(control_counts)

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'N sites', 'Value': sites_counts['Total']})
    results.append({'Metric': 'N control lysines', 'Value': control_counts['Total']})

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': '--- Overall category comparison ---', 'Value': ''})
    results.append({'Metric': 'Category', 'Value': 'Sites | Control'})

    for cat in ['Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic',
                'Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic']:
        s_count = sites_counts[cat]
        c_count = control_counts[cat]
        s_pct = 100 * s_count / sites_counts['Total'] if sites_counts['Total'] > 0 else 0
        c_pct = 100 * c_count / control_counts['Total'] if control_counts['Total'] > 0 else 0
        results.append({'Metric': cat, 'Value': f'{s_count} ({s_pct:.1f}%) | {c_count} ({c_pct:.1f}%)'})

    # Chi-square comparing overall distributions
    results.append({'Metric': '', 'Value': ''})
    sites_all = [sites_counts[c] for c in ['Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic',
                                            'Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic']]
    control_all = [control_counts[c] for c in ['Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic',
                                                'Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic']]

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
        ('Flexible', 'Flexible', 'flexible'),
        ('Consensus', None, None),  # Special handling
        ('Any_acidic', 'Any_acidic_within_threshold', 'any_acidic_within_threshold'),
        ('Exposed_acidic', 'Exposed_acidic_within_threshold', 'exposed_acidic_within_threshold'),
    ]

    for name, s_col_hint, c_col_hint in features:
        try:
            if name == 'Consensus':
                s_fwd = 'Forward_consensus' if 'Forward_consensus' in sites_df.columns else 'forward_consensus'
                s_inv = 'Inverse_consensus' if 'Inverse_consensus' in sites_df.columns else 'inverse_consensus'
                c_fwd = 'Forward_consensus' if 'Forward_consensus' in control_df.columns else 'forward_consensus'
                c_inv = 'Inverse_consensus' if 'Inverse_consensus' in control_df.columns else 'inverse_consensus'
                s_yes = ((sites_df[s_fwd] == 'Yes') | (sites_df[s_inv] == 'Yes')).sum()
                c_yes = ((control_df[c_fwd] == 'Yes') | (control_df[c_inv] == 'Yes')).sum()
            else:
                s_col = s_col_hint if s_col_hint in sites_df.columns else c_col_hint
                c_col = s_col_hint if s_col_hint in control_df.columns else c_col_hint
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
    """Analyze which features predict high score."""
    results = []
    results.append({'Metric': '=== PREDICTORS OF HIGH SCORE ===', 'Value': ''})
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

    # Define predictors - handle missing columns gracefully
    def get_mask(col, value='Yes'):
        if col in valid.columns:
            return valid[col] == value
        return pd.Series([False] * len(valid), index=valid.index)

    predictors = [
        ('Flexible', get_mask('Flexible')),
        ('Structured', get_mask('Structured')),
        ('Forward_consensus', get_mask('Forward_consensus')),
        ('Inverse_consensus', get_mask('Inverse_consensus')),
        ('Any_consensus', get_mask('Forward_consensus') | get_mask('Inverse_consensus')),
        ('Any_acidic_within_threshold', get_mask('Any_acidic_within_threshold')),
        ('Exposed_acidic_within_threshold', get_mask('Exposed_acidic_within_threshold')),
        ('Acidic_in_pm2', get_mask('Acidic_in_pm2')),
        ('Exposed_acidic_in_pm2', get_mask('Exposed_acidic_in_pm2')),
        ('Any_hydrophobic_within_threshold', get_mask('Any_hydrophobic_within_threshold')),
        ('Exposed_hydrophobic_within_threshold', get_mask('Exposed_hydrophobic_within_threshold')),
    ]

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

    # Multiple acidics as predictor (count from Acidics_within_threshold list)
    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': '--- N_acidic_within_threshold as predictor ---', 'Value': ''})

    acidics_col = 'Acidics_within_threshold'
    if acidics_col in valid.columns:
        # Count acidics from comma-separated list
        def count_acidics(val):
            if pd.isna(val) or val == '':
                return 0
            return len(str(val).split(','))

        valid = valid.copy()
        valid['n_acidic'] = valid[acidics_col].apply(count_acidics)

        # Correlation with score
        rho, p = stats.spearmanr(valid['n_acidic'], valid[score_col])
        results.append({'Metric': 'Spearman: N_acidic vs Score', 'Value': f'ρ={rho:.3f}, p={p:.2e}'})

        # Point-biserial with high_score
        r, p = stats.pointbiserialr(valid['high_score'], valid['n_acidic'])
        results.append({'Metric': 'Point-biserial: N_acidic vs high_score', 'Value': f'r={r:.3f}, p={p:.2e}'})

        # Mean N_acidic by score group
        results.append({'Metric': '', 'Value': ''})
        high_mean = valid[valid['high_score'] == 1]['n_acidic'].mean()
        low_mean = valid[valid['high_score'] == 0]['n_acidic'].mean()
        t_stat, t_p = stats.ttest_ind(
            valid[valid['high_score'] == 1]['n_acidic'],
            valid[valid['high_score'] == 0]['n_acidic']
        )
        results.append({'Metric': 'Mean N_acidic (high score)', 'Value': f'{high_mean:.2f}'})
        results.append({'Metric': 'Mean N_acidic (low score)', 'Value': f'{low_mean:.2f}'})
        results.append({'Metric': 't-test', 'Value': f't={t_stat:.2f}, p={t_p:.2e}'})

        # Test 2+ acidics
        results.append({'Metric': '', 'Value': ''})
        has_2plus = (valid['n_acidic'] >= 2).astype(int)
        tp = ((has_2plus == 1) & (valid['high_score'] == 1)).sum()
        fp = ((has_2plus == 1) & (valid['high_score'] == 0)).sum()
        fn = ((has_2plus == 0) & (valid['high_score'] == 1)).sum()
        tn = ((has_2plus == 0) & (valid['high_score'] == 0)).sum()

        if fp > 0 and fn > 0:
            odds = (tp * tn) / (fp * fn)
            _, p = stats.fisher_exact([[tp, fp], [fn, tn]])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            results.append({'Metric': 'Has 2+ acidics within threshold', 'Value': f'OR={odds:.2f}, p={p:.2e} {sig}'})

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

        # Define predictors (excluding Flexible/Structured since we've already filtered)
        def get_mask(col, value='Yes'):
            if col in subset.columns:
                return subset[col] == value
            return pd.Series([False] * len(subset), index=subset.index)

        predictors = [
            ('Forward_consensus', get_mask('Forward_consensus')),
            ('Inverse_consensus', get_mask('Inverse_consensus')),
            ('Any_consensus', get_mask('Forward_consensus') | get_mask('Inverse_consensus')),
            ('Any_acidic_within_threshold', get_mask('Any_acidic_within_threshold')),
            ('Exposed_acidic_within_threshold', get_mask('Exposed_acidic_within_threshold')),
            ('Acidic_in_pm2', get_mask('Acidic_in_pm2')),
            ('Exposed_acidic_in_pm2', get_mask('Exposed_acidic_in_pm2')),
            ('Any_hydrophobic_within_threshold', get_mask('Any_hydrophobic_within_threshold')),
            ('Exposed_hydrophobic_within_threshold', get_mask('Exposed_hydrophobic_within_threshold')),
        ]

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

        # N_acidic analysis (count from Acidics_within_threshold list)
        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': '--- N_acidic_within_threshold as predictor ---', 'Value': ''})

        acidics_col = 'Acidics_within_threshold'
        if acidics_col in subset.columns:
            # Count acidics from comma-separated list
            def count_acidics(val):
                if pd.isna(val) or val == '':
                    return 0
                return len(str(val).split(','))

            subset = subset.copy()
            subset['n_acidic'] = subset[acidics_col].apply(count_acidics)

            # Correlation with score
            if len(subset) >= 10:
                rho, p = stats.spearmanr(subset['n_acidic'], subset[score_col])
                results.append({'Metric': 'Spearman: N_acidic vs Score', 'Value': f'ρ={rho:.3f}, p={p:.2e}'})

                # Point-biserial with high_score
                r, p = stats.pointbiserialr(subset['high_score'], subset['n_acidic'])
                results.append({'Metric': 'Point-biserial: N_acidic vs high_score', 'Value': f'r={r:.3f}, p={p:.2e}'})

                # Mean N_acidic by score group
                results.append({'Metric': '', 'Value': ''})
                high_mean = subset[subset['high_score'] == 1]['n_acidic'].mean()
                low_mean = subset[subset['high_score'] == 0]['n_acidic'].mean()

                if n_high >= 2 and n_low >= 2:
                    t_stat, t_p = stats.ttest_ind(
                        subset[subset['high_score'] == 1]['n_acidic'],
                        subset[subset['high_score'] == 0]['n_acidic']
                    )
                    results.append({'Metric': 'Mean N_acidic (high score)', 'Value': f'{high_mean:.2f}'})
                    results.append({'Metric': 'Mean N_acidic (low score)', 'Value': f'{low_mean:.2f}'})
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


def analyze_acidic_distance_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze distribution of acidic distances (in residues) from lysine site.

    Calculates mean, median, std dev separately for:
    - Flexible sites
    - Structured sites (overall and broken down by helix/coil/strand)
    - Score groups (VeryHigh, High, Medium, Low)
    """
    results = []
    results.append({'Metric': '=== ACIDIC DISTANCE DISTRIBUTION (residues from Lysine) ===', 'Value': ''})
    results.append({'Metric': 'Note: Distance in residue positions (not Angstroms)', 'Value': ''})
    results.append({'Metric': '', 'Value': ''})

    # Column to get relative positions
    positions_col = 'Distances_positions'
    if positions_col not in df.columns:
        results.append({'Metric': 'Error', 'Value': 'Distances_positions column not found'})
        return pd.DataFrame(results)

    def extract_distances(val):
        """Extract numeric distances from position string like '-2(exp),+5,-10'."""
        if pd.isna(val) or val == '':
            return []
        distances = []
        for part in str(val).split(','):
            # Remove (exp) marker and convert to int
            part_clean = part.replace('(exp)', '').strip()
            try:
                d = int(part_clean)
                distances.append(abs(d))  # Use absolute distance
            except:
                pass
        return distances

    def calc_stats(distances_list):
        """Calculate mean, median, std for a list of distances."""
        if not distances_list:
            return None, None, None
        arr = np.array(distances_list)
        return np.mean(arr), np.median(arr), np.std(arr)

    def analyze_subset(subset_df, label):
        """Analyze a subset and return results list."""
        res = []
        distances = []
        for val in subset_df[positions_col].dropna():
            distances.extend(extract_distances(val))

        if distances:
            mean, median, std = calc_stats(distances)
            res.append({'Metric': f'--- {label} ---', 'Value': ''})
            res.append({'Metric': 'N sites', 'Value': len(subset_df)})
            res.append({'Metric': 'N acidics', 'Value': len(distances)})
            res.append({'Metric': 'Mean distance (residues)', 'Value': f'{mean:.2f}'})
            res.append({'Metric': 'Median distance', 'Value': f'{median:.1f}'})
            res.append({'Metric': 'Std dev', 'Value': f'{std:.2f}'})
            res.append({'Metric': '', 'Value': ''})
        return res

    # Overall stats
    results.extend(analyze_subset(df, 'ALL SITES'))

    # Flexible sites
    flex_df = df[df['Flexible'] == 'Yes']
    results.extend(analyze_subset(flex_df, 'FLEXIBLE SITES'))

    # Structured sites - overall
    struct_df = df[df['Structured'] == 'Yes']
    results.extend(analyze_subset(struct_df, 'STRUCTURED SITES (overall)'))

    # Structured sites by secondary structure
    ss_col = 'Secondary_structure'
    if ss_col in df.columns:
        for ss_type in ['helix', 'strand', 'coil']:
            ss_df = struct_df[struct_df[ss_col] == ss_type]
            results.extend(analyze_subset(ss_df, f'STRUCTURED - {ss_type.upper()}'))

    # By score groups
    score_col = 'Score (SUMO site)'
    if score_col in df.columns:
        df_with_score = df.copy()
        df_with_score[score_col] = pd.to_numeric(df_with_score[score_col], errors='coerce')

        results.append({'Metric': '========================================', 'Value': ''})
        results.append({'Metric': 'BY SCORE GROUPS', 'Value': ''})
        results.append({'Metric': '========================================', 'Value': ''})
        results.append({'Metric': '', 'Value': ''})

        score_groups = {
            'VeryHigh': df_with_score[df_with_score[score_col] >= SCORE_VERY_HIGH_THRESHOLD],
            'High': df_with_score[(df_with_score[score_col] >= SCORE_HIGH_THRESHOLD) & (df_with_score[score_col] < SCORE_VERY_HIGH_THRESHOLD)],
            'Medium': df_with_score[(df_with_score[score_col] >= SCORE_MEDIUM_THRESHOLD) & (df_with_score[score_col] < SCORE_HIGH_THRESHOLD)],
            'Low': df_with_score[df_with_score[score_col] < SCORE_MEDIUM_THRESHOLD]
        }

        for group_name, group_df in score_groups.items():
            if len(group_df) > 0:
                results.extend(analyze_subset(group_df, f'SCORE: {group_name}'))

    return pd.DataFrame(results)


def analyze_acidic_distance_angstroms(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze distribution of acidic distances (in Angstroms) from lysine site.

    Calculates mean, median, std dev separately for:
    - Flexible sites
    - Structured sites (overall and broken down by helix/coil/strand)
    - Score groups (VeryHigh, High, Medium, Low)
    """
    results = []
    results.append({'Metric': '=== ACIDIC DISTANCE DISTRIBUTION (Angstroms) ===', 'Value': ''})
    results.append({'Metric': 'Note: Distance in Angstroms (Cα-Cα)', 'Value': ''})
    results.append({'Metric': '', 'Value': ''})

    # Column to get Angstrom distances
    angstrom_col = 'Distances_Angstroms'
    if angstrom_col not in df.columns:
        results.append({'Metric': 'Error', 'Value': 'Distances_Angstroms column not found'})
        return pd.DataFrame(results)

    def extract_distances(val):
        """Extract numeric distances from angstrom string like '5.2(exp),6.8,7.1'."""
        if pd.isna(val) or val == '':
            return []
        distances = []
        for part in str(val).split(','):
            # Remove (exp) marker and convert to float
            part_clean = part.replace('(exp)', '').strip()
            try:
                d = float(part_clean)
                distances.append(d)
            except:
                pass
        return distances

    def calc_stats(distances_list):
        """Calculate mean, median, std for a list of distances."""
        if not distances_list:
            return None, None, None
        arr = np.array(distances_list)
        return np.mean(arr), np.median(arr), np.std(arr)

    def analyze_subset(subset_df, label):
        """Analyze a subset and return results list."""
        res = []
        distances = []
        for val in subset_df[angstrom_col].dropna():
            distances.extend(extract_distances(val))

        if distances:
            mean, median, std = calc_stats(distances)
            res.append({'Metric': f'--- {label} ---', 'Value': ''})
            res.append({'Metric': 'N sites', 'Value': len(subset_df)})
            res.append({'Metric': 'N acidics', 'Value': len(distances)})
            res.append({'Metric': 'Mean distance (Å)', 'Value': f'{mean:.2f}'})
            res.append({'Metric': 'Median distance (Å)', 'Value': f'{median:.2f}'})
            res.append({'Metric': 'Std dev', 'Value': f'{std:.2f}'})
            res.append({'Metric': '', 'Value': ''})
        return res

    # Overall stats
    results.extend(analyze_subset(df, 'ALL SITES'))

    # Flexible sites
    flex_df = df[df['Flexible'] == 'Yes']
    results.extend(analyze_subset(flex_df, 'FLEXIBLE SITES'))

    # Structured sites - overall
    struct_df = df[df['Structured'] == 'Yes']
    results.extend(analyze_subset(struct_df, 'STRUCTURED SITES (overall)'))

    # Structured sites by secondary structure
    ss_col = 'Secondary_structure'
    if ss_col in df.columns:
        for ss_type in ['helix', 'strand', 'coil']:
            ss_df = struct_df[struct_df[ss_col] == ss_type]
            results.extend(analyze_subset(ss_df, f'STRUCTURED - {ss_type.upper()}'))

    # By score groups
    score_col = 'Score (SUMO site)'
    if score_col in df.columns:
        df_with_score = df.copy()
        df_with_score[score_col] = pd.to_numeric(df_with_score[score_col], errors='coerce')

        results.append({'Metric': '========================================', 'Value': ''})
        results.append({'Metric': 'BY SCORE GROUPS', 'Value': ''})
        results.append({'Metric': '========================================', 'Value': ''})
        results.append({'Metric': '', 'Value': ''})

        score_groups = {
            'VeryHigh': df_with_score[df_with_score[score_col] >= SCORE_VERY_HIGH_THRESHOLD],
            'High': df_with_score[(df_with_score[score_col] >= SCORE_HIGH_THRESHOLD) & (df_with_score[score_col] < SCORE_VERY_HIGH_THRESHOLD)],
            'Medium': df_with_score[(df_with_score[score_col] >= SCORE_MEDIUM_THRESHOLD) & (df_with_score[score_col] < SCORE_HIGH_THRESHOLD)],
            'Low': df_with_score[df_with_score[score_col] < SCORE_MEDIUM_THRESHOLD]
        }

        for group_name, group_df in score_groups.items():
            if len(group_df) > 0:
                results.extend(analyze_subset(group_df, f'SCORE: {group_name}'))

    return pd.DataFrame(results)


def analyze_hydrophobic_features(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze hydrophobic features - percentages with hydrophobic within threshold.

    Shows statistics for:
    - All sites
    - Flexible sites
    - Structured sites
    """
    results = []
    results.append({'Metric': '=== HYDROPHOBIC ANALYSIS ===', 'Value': ''})
    results.append({'Metric': f'Hydrophobic distance threshold: {HYDROPHOBIC_DISTANCE_MIN} - {HYDROPHOBIC_DISTANCE_MAX} Å', 'Value': ''})
    results.append({'Metric': '', 'Value': ''})

    hp_col = 'Any_hydrophobic_within_threshold'
    hp_exp_col = 'Exposed_hydrophobic_within_threshold'

    if hp_col not in df.columns:
        results.append({'Metric': 'Error', 'Value': 'Hydrophobic columns not found'})
        return pd.DataFrame(results)

    valid_df = df[df['pLDDT_site'].notna()]

    def calc_hp_stats(subset, label):
        n_total = len(subset)
        if n_total == 0:
            return []
        n_hp = (subset[hp_col] == 'Yes').sum() if hp_col in subset.columns else 0
        n_hp_exp = (subset[hp_exp_col] == 'Yes').sum() if hp_exp_col in subset.columns else 0

        return [
            {'Metric': f'--- {label} ---', 'Value': ''},
            {'Metric': 'Total sites', 'Value': n_total},
            {'Metric': 'With hydrophobic within threshold', 'Value': f'{n_hp} ({100*n_hp/n_total:.1f}%)'},
            {'Metric': 'With exposed hydrophobic within threshold', 'Value': f'{n_hp_exp} ({100*n_hp_exp/n_total:.1f}%)'},
            {'Metric': '', 'Value': ''}
        ]

    results.extend(calc_hp_stats(valid_df, 'ALL SITES'))

    flex_df = valid_df[valid_df['Flexible'] == 'Yes']
    results.extend(calc_hp_stats(flex_df, 'FLEXIBLE SITES'))

    struct_df = valid_df[valid_df['Structured'] == 'Yes']
    results.extend(calc_hp_stats(struct_df, 'STRUCTURED SITES'))

    # Fisher test: Flexible vs Structured for hydrophobic
    if len(flex_df) > 0 and len(struct_df) > 0:
        results.append({'Metric': '--- STATISTICAL COMPARISON ---', 'Value': ''})

        flex_hp = (flex_df[hp_col] == 'Yes').sum() if hp_col in flex_df.columns else 0
        flex_no_hp = len(flex_df) - flex_hp
        struct_hp = (struct_df[hp_col] == 'Yes').sum() if hp_col in struct_df.columns else 0
        struct_no_hp = len(struct_df) - struct_hp

        try:
            odds, p = stats.fisher_exact([[flex_hp, flex_no_hp], [struct_hp, struct_no_hp]])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            flex_rate = flex_hp / len(flex_df) if len(flex_df) > 0 else 0
            struct_rate = struct_hp / len(struct_df) if len(struct_df) > 0 else 0
            results.append({
                'Metric': 'Any_hydrophobic: Flex vs Struct',
                'Value': f'Flex={flex_rate:.1%}, Struct={struct_rate:.1%}, OR={odds:.2f}, p={p:.2e} {sig}'
            })
        except Exception as e:
            results.append({'Metric': 'Fisher test error', 'Value': str(e)})

        # Same for exposed hydrophobic
        flex_hp_exp = (flex_df[hp_exp_col] == 'Yes').sum() if hp_exp_col in flex_df.columns else 0
        flex_no_hp_exp = len(flex_df) - flex_hp_exp
        struct_hp_exp = (struct_df[hp_exp_col] == 'Yes').sum() if hp_exp_col in struct_df.columns else 0
        struct_no_hp_exp = len(struct_df) - struct_hp_exp

        try:
            odds, p = stats.fisher_exact([[flex_hp_exp, flex_no_hp_exp], [struct_hp_exp, struct_no_hp_exp]])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            flex_rate = flex_hp_exp / len(flex_df) if len(flex_df) > 0 else 0
            struct_rate = struct_hp_exp / len(struct_df) if len(struct_df) > 0 else 0
            results.append({
                'Metric': 'Exposed_hydrophobic: Flex vs Struct',
                'Value': f'Flex={flex_rate:.1%}, Struct={struct_rate:.1%}, OR={odds:.2f}, p={p:.2e} {sig}'
            })
        except Exception as e:
            results.append({'Metric': 'Fisher test error', 'Value': str(e)})

    return pd.DataFrame(results)


def perform_global_analysis(df: pd.DataFrame, control_df: pd.DataFrame = None) -> pd.DataFrame:
    """Perform global analysis on all sites."""
    results = []

    # Header
    results.append({'Metric': '=== GLOBAL ANALYSIS ===', 'Value': ''})
    results.append({'Metric': '', 'Value': ''})

    # Category counts for sites
    counts = count_categories(df)
    rates = calculate_rates(counts)

    results.append({'Metric': '--- SUMO Sites Category Counts ---', 'Value': ''})
    for cat in ['Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic']:
        rate = rates.get(f'{cat}_rate', 0) or 0
        results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%} of Flexible)'})

    results.append({'Metric': '', 'Value': ''})
    for cat in ['Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic']:
        rate = rates.get(f'{cat}_rate', 0) or 0
        results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%} of Structured)'})

    results.append({'Metric': '', 'Value': ''})
    results.append({'Metric': 'Total Flexible', 'Value': counts['Total_Flexible']})
    results.append({'Metric': 'Total Structured', 'Value': counts['Total_Structured']})
    results.append({'Metric': 'Total', 'Value': counts['Total']})

    # Flex vs Struct comparison
    results.append({'Metric': '', 'Value': ''})
    results.extend(compare_flex_vs_struct(df, "SITES:"))

    # Control analysis
    if control_df is not None and len(control_df) > 0:
        # Check if control has category data
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
        for cat in ['Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic']:
            rate = control_rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{control_counts[cat]} ({rate:.1%} of Flexible)'})

        results.append({'Metric': '', 'Value': ''})
        for cat in ['Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic']:
            rate = control_rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{control_counts[cat]} ({rate:.1%} of Structured)'})

        results.append({'Metric': '', 'Value': ''})
        results.append({'Metric': 'Total Flexible', 'Value': control_counts['Total_Flexible']})
        results.append({'Metric': 'Total Structured', 'Value': control_counts['Total_Structured']})
        results.append({'Metric': 'Total', 'Value': control_counts['Total']})

        # Flex vs Struct for control
        results.append({'Metric': '', 'Value': ''})
        results.extend(compare_flex_vs_struct(control_df, "CONTROL:"))

        # Sites vs Control comparison
        results.append({'Metric': '', 'Value': ''})
        results.extend(compare_sites_vs_control(df, control_df))

    return pd.DataFrame(results)


def perform_stratified_analysis(df: pd.DataFrame, score_col: str = 'Score (SUMO site)') -> dict:
    """Perform analysis stratified by score groups."""
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
        for cat in ['Flexible_consensus', 'Flexible_exposed_acidic', 'Flexible_buried_acidic', 'Flexible_no_acidic']:
            rate = rates.get(f'{cat}_rate', 0) or 0
            results.append({'Metric': cat, 'Value': f'{counts[cat]} ({rate:.1%})'})

        results.append({'Metric': '', 'Value': ''})
        for cat in ['Structured_consensus', 'Structured_exposed_acidic', 'Structured_buried_acidic', 'Structured_no_acidic']:
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
    """Load core output and perform all analyses."""
    print(f"\n{'=' * 60}\nSUMO Site Analyzer - Statistical Analysis\n{'=' * 60}\n")
    print(f"Parameters:")
    print(f"  pLDDT threshold: {PLDDT_THRESHOLD}")
    print(f"  Acidic distance threshold (Cα-Cα): {DISTANCE_THRESHOLD_MIN} - {DISTANCE_THRESHOLD_MAX} Å")
    print(f"  Hydrophobic distance threshold (Cα-Cα): {HYDROPHOBIC_DISTANCE_MIN} - {HYDROPHOBIC_DISTANCE_MAX} Å")
    print(f"  Exposure: {EXPOSURE_DISTANCE_THRESHOLD} Å radius, <{EXPOSURE_NEIGHBOR_THRESHOLD} neighbors = exposed")
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
        output_file = os.path.splitext(input_file)[0] + '_analysis.xlsx'

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

        # Acidic distance distribution analysis (in residue positions)
        acidic_dist_df = analyze_acidic_distance_distribution(valid_df)
        acidic_dist_df.to_excel(writer, sheet_name='Acidic_Dist_Residues', index=False)

        # Acidic distance distribution analysis (in Angstroms)
        acidic_angstrom_df = analyze_acidic_distance_angstroms(valid_df)
        acidic_angstrom_df.to_excel(writer, sheet_name='Acidic_Dist_Angstroms', index=False)

        # Hydrophobic analysis
        hp_analysis_df = analyze_hydrophobic_features(valid_df)
        hp_analysis_df.to_excel(writer, sheet_name='Hydrophobic_Analysis', index=False)

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
        print("\nConfigurable parameters (must match sumo_site_core.py):")
        print(f"  PLDDT_THRESHOLD = {PLDDT_THRESHOLD}")
        print(f"  DISTANCE_THRESHOLD_MIN = {DISTANCE_THRESHOLD_MIN} Å (Cα-Cα, acidic)")
        print(f"  DISTANCE_THRESHOLD_MAX = {DISTANCE_THRESHOLD_MAX} Å (Cα-Cα, acidic)")
        print(f"  HYDROPHOBIC_DISTANCE_MIN = {HYDROPHOBIC_DISTANCE_MIN} Å (Cα-Cα)")
        print(f"  HYDROPHOBIC_DISTANCE_MAX = {HYDROPHOBIC_DISTANCE_MAX} Å (Cα-Cα)")
        print(f"  EXPOSURE_DISTANCE_THRESHOLD = {EXPOSURE_DISTANCE_THRESHOLD} Å")
        print(f"  EXPOSURE_NEIGHBOR_THRESHOLD = {EXPOSURE_NEIGHBOR_THRESHOLD}")
        print(f"  SCORE_VERY_HIGH_THRESHOLD = {SCORE_VERY_HIGH_THRESHOLD}")
        print(f"  SCORE_HIGH_THRESHOLD = {SCORE_HIGH_THRESHOLD}")
        print(f"  SCORE_MEDIUM_THRESHOLD = {SCORE_MEDIUM_THRESHOLD}")
        sys.exit(1)
    analyze_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
