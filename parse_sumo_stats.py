#!/usr/bin/env python3
"""
SUMO Site Analyzer (Stats Extended)
===================================
Analyzes SUMO sites using AlphaFold models.
- Categorizes sites: Consensus, Flank_acidic, Space_acidic, Neither.
- Measures distances using NZ (Lys) to CG (Asp) / CD (Glu).
- Advanced statistical analyses: logistic regression, ROC, bootstrap CIs.
- Random lysine background comparison.

Usage:
    python parse_sumo_stats.py <input_excel_file> [output_file]
"""

import sys
import os
import re
import math
import random
import requests
import pandas as pd
import numpy as np
from scipy import stats
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

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

# Parallel processing settings
MAX_WORKERS = 10

# Bootstrap settings
N_BOOTSTRAP = 1000

# =============================================================================

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Global session for connection pooling
_session = None
_session_lock = Lock()

# Store structures for random lysine analysis
_structure_cache = {}
_structure_cache_lock = Lock()


def get_session() -> requests.Session:
    global _session
    if _session is None:
        with _session_lock:
            if _session is None:
                _session = requests.Session()
                adapter = requests.adapters.HTTPAdapter(pool_connections=MAX_WORKERS, pool_maxsize=MAX_WORKERS * 2)
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


def parse_uniprot_id(protein_id: str) -> tuple[str, int | None]:
    match = re.search(r'([A-Z][0-9][A-Z0-9]{3,9}[0-9])(?:-(\d+))?', protein_id)
    if match:
        return match.group(1), int(match.group(2)) if match.group(2) else None
    return protein_id, None


def safe_get(url: str, timeout: int = 30) -> requests.Response | None:
    try:
        r = get_session().get(url, timeout=timeout)
        return r if r.status_code == 200 else None
    except Exception:
        return None


def calc_distance(c1: tuple, c2: tuple) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


# =============================================================================
# SEQUENCE & MAPPING LOGIC
# =============================================================================

def get_unisave_sequence(uid: str, cache: dict, min_length: int = 0) -> str | None:
    key = f"unisave_{uid}" if min_length == 0 else f"unisave_{uid}_min{min_length}"
    if key in cache:
        return cache[key]

    versions_to_try = [1, 2, 3, 5, 10, 20, 30, 40, 50]
    best_seq = None

    for version in versions_to_try:
        r = safe_get(f"https://rest.uniprot.org/unisave/{uid}?format=fasta&versions={version}")
        if r and r.text.strip() and not r.text.startswith('<!'):
            lines = r.text.strip().split('\n')
            seq = ''.join(l for l in lines if not l.startswith('>'))
            if seq and len(seq) > 10:
                if min_length > 0:
                    if len(seq) >= min_length:
                        cache[key] = seq
                        return seq
                    if best_seq is None or len(seq) > len(best_seq):
                        best_seq = seq
                else:
                    cache[key] = seq
                    return seq

    if best_seq:
        cache[key] = best_seq
        return best_seq

    cache[key] = None
    return None


def get_uniprot_sequence(uid: str, isoform: int | None, cache: dict) -> str | None:
    key = f"upseq_{uid}_{isoform}"
    if key in cache:
        return cache[key]
    suffix = f"-{isoform}" if isoform else ""
    r = safe_get(f"https://rest.uniprot.org/uniprotkb/{uid}{suffix}.fasta")
    if r:
        lines = r.text.strip().split('\n')
        seq = ''.join(l for l in lines if not l.startswith('>'))
        cache[key] = seq if seq else None
    else:
        cache[key] = None
    return cache[key]


def get_uniprot_info(uid: str, cache: dict) -> dict | None:
    key = f"upinfo_{uid}"
    if key in cache:
        return cache[key]
    r = safe_get(f"https://rest.uniprot.org/uniprotkb/{uid}.json")
    if r:
        try:
            data = r.json()
            genes = [g['geneName']['value'] for g in data.get('genes', []) if 'geneName' in g]
            cache[key] = {
                'genes': genes,
                'sequence': data.get('sequence', {}).get('value', ''),
                'primary': data.get('primaryAccession', uid)
            }
        except Exception:
            cache[key] = None
    else:
        cache[key] = None
    return cache[key]


def get_mapped_accession(uid: str, cache: dict) -> str | None:
    key = f"mapped_{uid}"
    if key in cache:
        return cache[key]
    info = get_uniprot_info(uid, cache)
    if info and info['primary'] and info['primary'] != uid:
        cache[key] = info['primary']
        return cache[key]
    r = safe_get(f"https://rest.uniprot.org/uniprotkb/search?query=(sec_acc:{uid})&fields=accession&size=1")
    if r:
        try:
            results = r.json().get('results', [])
            if results:
                acc = results[0].get('primaryAccession')
                if acc:
                    cache[key] = acc
                    return acc
        except Exception:
            pass
    cache[key] = None
    return None


def map_position_by_residue(src_seq: str, tgt_seq: str, src_pos: int) -> int | None:
    if not src_seq or not tgt_seq or src_pos < 1 or src_pos > len(src_seq):
        return None
    target_aa = src_seq[src_pos - 1]
    for w in [20, 15, 10, 7, 5]:
        start, end = max(0, src_pos - 1 - w), min(len(src_seq), src_pos + w)
        query, pos_in_q = src_seq[start:end], src_pos - 1 - start
        best_pos, best_score = None, 0
        for i in range(max(0, len(tgt_seq) - len(query) + 1)):
            tgt_window = tgt_seq[i:i + len(query)]
            if len(tgt_window) < len(query):
                continue
            matches = sum(a == b for a, b in zip(query, tgt_window))
            tgt_idx = i + pos_in_q
            if tgt_idx < len(tgt_seq) and tgt_seq[tgt_idx] == target_aa:
                score = matches + 10
                if score > best_score:
                    best_score, best_pos = score, tgt_idx + 1
        if best_pos and best_score >= len(query) * 0.4 + 10:
            return best_pos
    return None


# =============================================================================
# ALPHAFOLD HANDLING
# =============================================================================

def fetch_alphafold(uid: str, fragment: int = 1) -> dict | None:
    r = safe_get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}")
    if not r:
        return None
    try:
        data = r.json()
        if not data:
            return None
        cif_url = next((e['cifUrl'] for e in data if f"-F{fragment}-" in e.get('cifUrl', '')), None)
        if not cif_url and fragment == 1 and data:
            cif_url = data[0].get('cifUrl')
        if not cif_url:
            return None
        cif_r = safe_get(cif_url)
        return parse_cif(cif_r.text) if cif_r else None
    except Exception:
        return None


def parse_cif(content: str) -> dict:
    """Parse mmCIF file to extract coordinates and pLDDT values."""
    plddt, seq, coords_ca, coords_nz, acidic_coords = {}, {}, {}, {}, {}
    lysine_positions = []
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
                    if aa == 'K':
                        lysine_positions.append(rn)
                elif atom_name == 'NZ' and aa3 == 'LYS':
                    coords_nz[rn] = coord
                elif atom_name == 'CG' and aa3 == 'ASP':
                    acidic_coords[rn] = coord
                elif atom_name == 'CD' and aa3 == 'GLU':
                    acidic_coords[rn] = coord
            except Exception:
                continue
        elif in_atom and line.startswith('#'):
            break

    seq_str = ''.join(seq.get(i, 'X') for i in range(1, max(seq.keys(), default=0) + 1))
    acidic_positions = sorted(acidic_coords.keys())

    return {
        'plddt': plddt,
        'sequence': seq,
        'sequence_string': seq_str,
        'coords_ca': coords_ca,
        'coords_nz': coords_nz,
        'acidic_coords': acidic_coords,
        'acidic_positions': acidic_positions,
        'lysine_positions': lysine_positions
    }


# =============================================================================
# RESOLUTION & ANALYSIS
# =============================================================================

def resolve_structure_and_position(uid: str, isoform: int | None, pos: int, cache: dict) -> tuple[dict | None, str, int]:
    cache_key = f"resolve_{uid}_{isoform}_{pos}"
    if cache_key in cache:
        return cache[cache_key]

    struct = fetch_alphafold(uid, 1)
    if struct and pos > max(struct['plddt'].keys(), default=0):
        for f in range(2, 21):
            next_struct = fetch_alphafold(uid, f)
            if not next_struct:
                break
            if pos in next_struct['plddt']:
                struct = next_struct
                break

    if struct:
        # Store structure for random lysine analysis
        with _structure_cache_lock:
            _structure_cache[uid] = struct

        final_pos = pos
        if isoform:
            iso_seq = get_uniprot_sequence(uid, isoform, cache)
            if iso_seq:
                mapped = map_position_by_residue(iso_seq, struct['sequence_string'], pos)
                if mapped:
                    final_pos = mapped
                else:
                    return (None, uid, pos)
        result = (struct, uid, final_pos)
        cache[cache_key] = result
        return result

    original_seq = get_unisave_sequence(uid, cache, min_length=pos) or get_uniprot_sequence(uid, isoform, cache)
    if not original_seq or pos > len(original_seq):
        return (None, uid, pos)

    target_acc = get_mapped_accession(uid, cache)
    if target_acc:
        struct = fetch_alphafold(target_acc, 1)
        if struct:
            with _structure_cache_lock:
                _structure_cache[target_acc] = struct
            mapped_pos = map_position_by_residue(original_seq, struct['sequence_string'], pos)
            if mapped_pos:
                return (struct, target_acc, mapped_pos)

    return (None, uid, pos)


def analyze_site(struct: dict, pos: int) -> dict:
    """Analyze a SUMO site with extended metrics including acidic count."""
    res = {
        'plddt_site': None, 'plddt_11avg': None,
        'aa_m2': None, 'aa_m1': None, 'aa_p1': None, 'aa_p2': None,
        'dist_acidic_any': None, 'closest_acidic_any': None,
        'dist_acidic_space': None, 'closest_acidic_space': None,
        'n_acidic_within_threshold': 0,  # NEW: count of acidic residues within threshold
        'forward_consensus': None, 'inverse_consensus': None,
        'flexible': None, 'structured': None,
        'acidic_in_space': None, 'acidic_in_flank': None, 'category': None
    }

    if not struct or pos not in struct['plddt']:
        return res

    plddt = struct['plddt']
    seq = struct['sequence']
    coords_nz = struct['coords_nz']
    acidic_coords = struct['acidic_coords']
    acidic_positions = struct['acidic_positions']

    res['plddt_site'] = plddt.get(pos)
    window_vals = [plddt[i] for i in range(pos - 5, pos + 6) if i in plddt]
    res['plddt_11avg'] = sum(window_vals) / len(window_vals) if window_vals else None

    res['aa_m2'] = seq.get(pos - 2)
    res['aa_m1'] = seq.get(pos - 1)
    res['aa_p1'] = seq.get(pos + 1)
    res['aa_p2'] = seq.get(pos + 2)

    # Calculate distances using NZ of Lys to CG/CD of acidic residues
    if pos in coords_nz and acidic_positions:
        lys_nz_coord = coords_nz[pos]
        min_any, cls_any = float('inf'), None
        min_space, cls_space = float('inf'), None
        n_within = 0

        for ac in acidic_positions:
            if ac not in acidic_coords:
                continue
            d = calc_distance(lys_nz_coord, acidic_coords[ac])

            # Count acidic within threshold (excluding -2/+2)
            rel_pos = ac - pos
            if rel_pos != -2 and rel_pos != 2 and d <= DISTANCE_THRESHOLD:
                n_within += 1

            if d < min_any:
                min_any, cls_any = d, ac

            if rel_pos != -2 and rel_pos != 2:
                if d < min_space:
                    min_space, cls_space = d, ac

        res['n_acidic_within_threshold'] = n_within

        if cls_any:
            res['dist_acidic_any'] = round(min_any, 2)
            res['closest_acidic_any'] = f"{seq.get(cls_any, '?')}{cls_any}"
        if cls_space:
            res['dist_acidic_space'] = round(min_space, 2)
            res['closest_acidic_space'] = f"{seq.get(cls_space, '?')}{cls_space}"

    # Consensus motifs
    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    res['inverse_consensus'] = 'Yes' if res['aa_p1'] in HYDROPHOBIC_RESIDUES and res['aa_m2'] in ACIDIC_RESIDUES else None

    # Flexibility classification
    if res['plddt_11avg'] is not None:
        if res['plddt_11avg'] < PLDDT_THRESHOLD:
            res['flexible'] = 'Yes'
        else:
            res['structured'] = 'Yes'

    # Acidic in space
    if res['dist_acidic_space'] is not None and res['dist_acidic_space'] <= DISTANCE_THRESHOLD:
        res['acidic_in_space'] = 'Yes'

    # Acidic in flank
    if res['aa_m2'] in ACIDIC_RESIDUES or res['aa_p2'] in ACIDIC_RESIDUES:
        res['acidic_in_flank'] = 'Yes'

    # Category assignment
    is_consensus = (res['forward_consensus'] == 'Yes' or res['inverse_consensus'] == 'Yes')
    is_flank_acidic = (res['acidic_in_flank'] == 'Yes')
    is_space_acidic = (res['acidic_in_space'] == 'Yes')

    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"

    if is_consensus:
        suffix = "_consensus"
    elif is_flank_acidic:
        suffix = "_flank_acidic"
    elif is_space_acidic:
        suffix = "_space_acidic"
    else:
        suffix = "_neither"

    res['category'] = prefix + suffix
    return res


# =============================================================================
# ADVANCED STATISTICAL ANALYSES
# =============================================================================

def bootstrap_rate_ci(successes: int, total: int, n_bootstrap: int = N_BOOTSTRAP, ci: float = 0.95) -> tuple:
    """Calculate bootstrap confidence interval for a rate."""
    if total == 0:
        return (0.0, 0.0, 0.0)
    rate = successes / total
    if total < 5:
        return (rate, rate, rate)

    # Generate bootstrap samples
    samples = np.random.binomial(total, rate, n_bootstrap) / total
    alpha = (1 - ci) / 2
    lower = np.percentile(samples, alpha * 100)
    upper = np.percentile(samples, (1 - alpha) * 100)
    return (rate, lower, upper)


def perform_logistic_regression(df: pd.DataFrame) -> pd.DataFrame:
    """
    Logistic regression: predict high score from features.
    Uses simple approach without sklearn dependency.
    """
    data_list = []
    data_list.append({'Metric': '=== LOGISTIC REGRESSION ANALYSIS ===', 'Value': ''})
    data_list.append({'Metric': 'Predicting: Is site high-scoring (>=400)?', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    score_col = 'Score (SUMO site)'
    valid = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    if len(valid) < 50:
        data_list.append({'Metric': 'Error', 'Value': 'Too few samples for regression'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    # Create binary outcome and features
    valid['high_score'] = (valid[score_col] >= SCORE_HIGH_THRESHOLD).astype(int)
    valid['is_flexible'] = (valid['Flexible'] == 'Yes').astype(int)
    valid['is_consensus'] = ((valid['Forward_consensus'] == 'Yes') | (valid['Inverse_consensus'] == 'Yes')).astype(int)
    valid['has_flank'] = (valid['Acidic_in_-2/+2'] == 'Yes').astype(int)
    valid['has_space'] = (valid['Acidic_in_space'] == 'Yes').astype(int)
    valid['has_any_acidic'] = ((valid['has_flank'] == 1) | (valid['has_space'] == 1)).astype(int)

    # Use distance as continuous variable
    valid['dist_space'] = pd.to_numeric(valid['Dist_acidic_space_A'], errors='coerce').fillna(50)

    features = ['is_flexible', 'is_consensus', 'has_flank', 'has_space']

    # Calculate univariate odds ratios
    data_list.append({'Metric': '--- Univariate Odds Ratios ---', 'Value': ''})

    for feat in features:
        # 2x2 table
        a = len(valid[(valid[feat] == 1) & (valid['high_score'] == 1)])
        b = len(valid[(valid[feat] == 1) & (valid['high_score'] == 0)])
        c = len(valid[(valid[feat] == 0) & (valid['high_score'] == 1)])
        d = len(valid[(valid[feat] == 0) & (valid['high_score'] == 0)])

        if b > 0 and c > 0 and d > 0:
            odds_ratio = (a * d) / (b * c) if b * c > 0 else float('inf')
            # Fisher exact test
            try:
                _, p = stats.fisher_exact([[a, b], [c, d]])
                sig = '*' if p < 0.05 else ''
                data_list.append({'Metric': f'{feat}', 'Value': f'OR={odds_ratio:.2f}, p={p:.2e} {sig}'})
            except:
                data_list.append({'Metric': f'{feat}', 'Value': f'OR={odds_ratio:.2f}'})

    data_list.append({'Metric': '', 'Value': ''})

    # Distance as continuous predictor (point-biserial correlation)
    data_list.append({'Metric': '--- Distance as Continuous Predictor ---', 'Value': ''})

    # Correlation between distance and high_score
    valid_dist = valid[valid['dist_space'] < 50]  # Exclude missing
    if len(valid_dist) > 20:
        rho, p = stats.pointbiserialr(valid_dist['high_score'], valid_dist['dist_space'])
        data_list.append({'Metric': 'Distance vs high_score correlation', 'Value': f'r={rho:.3f}, p={p:.2e}'})

        # Mean distance comparison
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

    # Summary statistics
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

    # Mann-Whitney U test (non-parametric)
    u_stat, u_p = stats.mannwhitneyu(flex, struct, alternative='two-sided')
    data_list.append({'Metric': 'Mann-Whitney U test', 'Value': f'U={u_stat:.0f}, p={u_p:.2e}'})

    # T-test
    t_stat, t_p = stats.ttest_ind(flex, struct)
    data_list.append({'Metric': 't-test', 'Value': f't={t_stat:.2f}, p={t_p:.2e}'})

    # Kolmogorov-Smirnov test
    ks_stat, ks_p = stats.ks_2samp(flex, struct)
    data_list.append({'Metric': 'KS test (distribution shape)', 'Value': f'D={ks_stat:.3f}, p={ks_p:.2e}'})

    # Distance percentiles
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

    # Correlations
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

    # Mean pLDDT by score category
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

    # For each predictor, calculate sensitivity, specificity, PPV, NPV
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

    # AUC approximation using Mann-Whitney U
    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- AUC estimates (distance as continuous) ---', 'Value': ''})

    valid['dist'] = pd.to_numeric(valid['Dist_acidic_space_A'], errors='coerce')
    valid_dist = valid[valid['dist'].notna()]

    if len(valid_dist) > 20:
        # AUC = P(random positive has lower distance than random negative)
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

        # Count acidic at this position
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

        # Fisher exact test
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

    # Rates to bootstrap
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

    # Bootstrap difference in rates
    if len(flex) > 10 and len(struct) > 10:
        # Cons+flank rate difference
        flex_cf = len(flex[flex['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic'])]) / len(flex)
        struct_cf = len(struct[struct['Category'].isin(['Structured_consensus', 'Structured_flank_acidic'])]) / len(struct)
        obs_diff = flex_cf - struct_cf

        # Bootstrap the difference
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

    # Distribution of counts
    data_list.append({'Metric': '--- Distribution of acidic counts ---', 'Value': ''})
    for n in range(5):
        count = len(valid[valid['n_acidic'] == n])
        pct = 100 * count / len(valid) if len(valid) > 0 else 0
        data_list.append({'Metric': f'{n} acidic residues', 'Value': f'{count} ({pct:.1f}%)'})

    count_3plus = len(valid[valid['n_acidic'] >= 3])
    pct_3plus = 100 * count_3plus / len(valid) if len(valid) > 0 else 0
    data_list.append({'Metric': '3+ acidic residues', 'Value': f'{count_3plus} ({pct_3plus:.1f}%)'})

    data_list.append({'Metric': '', 'Value': ''})

    # Score by acidic count
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

        # Correlation between count and score
        rho, p = stats.spearmanr(valid['n_acidic'], valid[score_col])
        data_list.append({'Metric': 'Spearman: n_acidic vs score', 'Value': f'rho={rho:.3f}, p={p:.2e}'})

    # Compare 0 vs 1+ acidic
    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- Comparing 0 vs 1+ acidic ---', 'Value': ''})

    has_0 = valid[valid['n_acidic'] == 0]
    has_1plus = valid[valid['n_acidic'] >= 1]

    if len(has_0) > 10 and len(has_1plus) > 10 and score_col in valid.columns:
        t_stat, p = stats.ttest_ind(has_0[score_col].dropna(), has_1plus[score_col].dropna())
        data_list.append({'Metric': 't-test (score)', 'Value': f't={t_stat:.2f}, p={p:.2e}'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


def perform_random_lysine_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Compare SUMO sites to random lysines from the same proteins."""
    data_list = []
    data_list.append({'Metric': '=== RANDOM LYSINE BACKGROUND COMPARISON ===', 'Value': ''})
    data_list.append({'Metric': 'Comparing SUMO sites to non-SUMO lysines from same proteins', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    # Get SUMO site positions per protein
    valid = df[df['pLDDT_site'].notna()].copy()
    sumo_sites = set()
    for _, row in valid.iterrows():
        uid = row.get('UniProt_ID_used')
        pos = row.get('Mapped_position') if pd.notna(row.get('Mapped_position')) else None
        if not pos:
            # Try to get original position
            pos_col = next((c for c in df.columns if 'position' in c.lower()), None)
            if pos_col:
                pos = row.get(pos_col)
        if uid and pos:
            sumo_sites.add((uid, int(pos)))

    data_list.append({'Metric': 'SUMO sites analyzed', 'Value': len(sumo_sites)})
    data_list.append({'Metric': 'Proteins with structures', 'Value': len(_structure_cache)})

    if len(_structure_cache) < 5:
        data_list.append({'Metric': 'Warning', 'Value': 'Too few structures for random lysine analysis'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    # Collect random lysines (non-SUMO sites)
    random_lys_results = []
    sumo_results = []

    for uid, struct in _structure_cache.items():
        lysine_positions = struct.get('lysine_positions', [])
        coords_nz = struct.get('coords_nz', {})
        acidic_coords = struct.get('acidic_coords', {})
        acidic_positions = struct.get('acidic_positions', [])
        seq = struct.get('sequence', {})
        plddt = struct.get('plddt', {})

        for lys_pos in lysine_positions:
            is_sumo = (uid, lys_pos) in sumo_sites

            # Calculate properties for this lysine
            if lys_pos not in coords_nz or lys_pos not in plddt:
                continue

            lys_nz = coords_nz[lys_pos]

            # Find closest acidic (not at -2/+2)
            min_dist = float('inf')
            has_acidic_in_space = False
            n_acidic_within = 0

            for ac in acidic_positions:
                if ac not in acidic_coords:
                    continue
                rel = ac - lys_pos
                if rel == -2 or rel == 2:
                    continue
                d = calc_distance(lys_nz, acidic_coords[ac])
                if d < min_dist:
                    min_dist = d
                if d <= DISTANCE_THRESHOLD:
                    n_acidic_within += 1
                    has_acidic_in_space = True

            # Check flanking acidic
            has_flank = seq.get(lys_pos - 2) in ACIDIC_RESIDUES or seq.get(lys_pos + 2) in ACIDIC_RESIDUES

            # Check consensus
            aa_m1 = seq.get(lys_pos - 1)
            aa_p1 = seq.get(lys_pos + 1)
            aa_m2 = seq.get(lys_pos - 2)
            aa_p2 = seq.get(lys_pos + 2)
            is_consensus = (aa_m1 in HYDROPHOBIC_RESIDUES and aa_p2 in ACIDIC_RESIDUES) or \
                           (aa_p1 in HYDROPHOBIC_RESIDUES and aa_m2 in ACIDIC_RESIDUES)

            # pLDDT
            window_vals = [plddt[i] for i in range(lys_pos - 5, lys_pos + 6) if i in plddt]
            plddt_avg = sum(window_vals) / len(window_vals) if window_vals else None
            is_flexible = plddt_avg is not None and plddt_avg < PLDDT_THRESHOLD

            result = {
                'is_sumo': is_sumo,
                'has_acidic_in_space': has_acidic_in_space,
                'has_flank': has_flank,
                'is_consensus': is_consensus,
                'n_acidic': n_acidic_within,
                'min_dist': min_dist if min_dist < float('inf') else None,
                'is_flexible': is_flexible,
                'plddt_avg': plddt_avg
            }

            if is_sumo:
                sumo_results.append(result)
            else:
                random_lys_results.append(result)

    data_list.append({'Metric': 'Random lysines collected', 'Value': len(random_lys_results)})
    data_list.append({'Metric': 'SUMO lysines matched', 'Value': len(sumo_results)})
    data_list.append({'Metric': '', 'Value': ''})

    if len(random_lys_results) < 50 or len(sumo_results) < 20:
        data_list.append({'Metric': 'Warning', 'Value': 'Not enough lysines for comparison'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    # Compare rates
    data_list.append({'Metric': '--- Rate comparisons ---', 'Value': ''})

    for prop, label in [
        ('is_consensus', 'Consensus'),
        ('has_flank', 'Flank acidic'),
        ('has_acidic_in_space', 'Space acidic'),
    ]:
        sumo_rate = sum(1 for r in sumo_results if r[prop]) / len(sumo_results)
        random_rate = sum(1 for r in random_lys_results if r[prop]) / len(random_lys_results)

        # Fisher test
        a = sum(1 for r in sumo_results if r[prop])
        b = len(sumo_results) - a
        c = sum(1 for r in random_lys_results if r[prop])
        d = len(random_lys_results) - c

        try:
            odds, p = stats.fisher_exact([[a, b], [c, d]])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            data_list.append({
                'Metric': label,
                'Value': f'SUMO={sumo_rate:.3f}, Random={random_rate:.3f}, OR={odds:.2f}, p={p:.2e} {sig}'
            })
        except:
            data_list.append({
                'Metric': label,
                'Value': f'SUMO={sumo_rate:.3f}, Random={random_rate:.3f}'
            })

    # Compare flexibility
    data_list.append({'Metric': '', 'Value': ''})
    sumo_flex_rate = sum(1 for r in sumo_results if r['is_flexible']) / len(sumo_results)
    random_flex_rate = sum(1 for r in random_lys_results if r['is_flexible']) / len(random_lys_results)
    data_list.append({
        'Metric': 'Flexibility rate',
        'Value': f'SUMO={sumo_flex_rate:.3f}, Random={random_flex_rate:.3f}'
    })

    # Compare mean distance
    sumo_dists = [r['min_dist'] for r in sumo_results if r['min_dist'] is not None]
    random_dists = [r['min_dist'] for r in random_lys_results if r['min_dist'] is not None]

    if len(sumo_dists) > 10 and len(random_dists) > 10:
        data_list.append({'Metric': '', 'Value': ''})
        data_list.append({'Metric': 'Mean distance to acidic', 'Value': f'SUMO={np.mean(sumo_dists):.1f}Å, Random={np.mean(random_dists):.1f}Å'})
        t_stat, p = stats.ttest_ind(sumo_dists, random_dists)
        data_list.append({'Metric': 't-test', 'Value': f't={t_stat:.2f}, p={p:.2e}'})

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


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

    def rate_str(count, total):
        rate = count / total if total > 0 else 0
        return f"{rate:.3f} ({count}/{total})"

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
    """Global trend analysis."""
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
# MAIN PROCESSING
# =============================================================================

def process_file(input_file: str, output_file: str = None):
    print(f"\n{'=' * 60}\nSUMO Site Analyzer (Extended)\n{'=' * 60}\n")
    df = pd.read_excel(input_file, sheet_name=0, header=1)
    df = df.iloc[:, :27]

    protein_col = next((c for c in df.columns if str(c).lower() == 'protein'), 'Protein')
    position_col = next((c for c in df.columns if str(c).lower() == 'position'), 'Position')

    cache = {}

    def process_single_site(idx, row):
        uid, iso = parse_uniprot_id(str(row[protein_col]).strip())
        try:
            pos = int(row[position_col])
        except Exception:
            return idx, {'used_id': uid, 'mapped_pos': None}

        struct, used_id, final_pos = resolve_structure_and_position(uid, iso, pos, cache)
        analysis = analyze_site(struct, final_pos)
        analysis.update({'used_id': used_id, 'mapped_pos': final_pos if final_pos != pos else None})
        return idx, analysis

    results = [None] * len(df)
    completed = 0
    total = len(df)

    print(f"Processing {total} sites with {MAX_WORKERS} workers...")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_single_site, idx, row): idx for idx, row in df.iterrows()}

        for future in as_completed(futures):
            idx, analysis = future.result()
            results[idx] = analysis
            completed += 1
            if completed % 100 == 0 or completed == total:
                logger.info(f"Processed {completed}/{total} ({100 * completed / total:.1f}%)")

    res_df = pd.DataFrame(results)
    cols_map = {
        'used_id': 'UniProt_ID_used',
        'mapped_pos': 'Mapped_position',
        'plddt_site': 'pLDDT_site',
        'plddt_11avg': 'pLDDT_11residue_avg',
        'aa_m2': 'AA_minus2',
        'aa_m1': 'AA_minus1',
        'aa_p1': 'AA_plus1',
        'aa_p2': 'AA_plus2',
        'dist_acidic_any': 'Dist_acidic_any_A',
        'closest_acidic_any': 'Closest_acidic_any',
        'dist_acidic_space': 'Dist_acidic_space_A',
        'closest_acidic_space': 'Closest_acidic_space',
        'n_acidic_within_threshold': 'N_acidic_in_space',
        'forward_consensus': 'Forward_consensus',
        'inverse_consensus': 'Inverse_consensus',
        'flexible': 'Flexible',
        'structured': 'Structured',
        'acidic_in_space': 'Acidic_in_space',
        'acidic_in_flank': 'Acidic_in_-2/+2',
        'category': 'Category'
    }
    res_df = res_df.rename(columns=cols_map)

    ordered_cols = [
        'UniProt_ID_used', 'Mapped_position', 'pLDDT_site', 'pLDDT_11residue_avg',
        'AA_minus2', 'AA_minus1', 'AA_plus1', 'AA_plus2',
        'Dist_acidic_any_A', 'Closest_acidic_any', 'Dist_acidic_space_A', 'Closest_acidic_space',
        'N_acidic_in_space', 'Forward_consensus', 'Inverse_consensus', 'Flexible', 'Structured',
        'Acidic_in_space', 'Acidic_in_-2/+2', 'Category'
    ]
    res_df = res_df[[c for c in ordered_cols if c in res_df.columns]]

    df = pd.concat([df, res_df], axis=1)

    total_sites = len(df)
    success_count = df['pLDDT_site'].notna().sum()
    failed_count = total_sites - success_count
    success_rate = (success_count / total_sites * 100) if total_sites > 0 else 0

    print(f"\n{'=' * 40}")
    print("DATA COMPLETENESS SUMMARY")
    print(f"{'=' * 40}")
    print(f"Total sites:          {total_sites}")
    print(f"Successfully resolved: {success_count} ({success_rate:.1f}%)")
    print(f"Failed to resolve:     {failed_count}")
    print(f"{'=' * 40}\n")

    summary_stats = {'total': total_sites, 'success': success_count, 'failed': failed_count}

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
        output_file = os.path.splitext(input_file)[0] + '_analyzed_stats.xlsx'

    print("Generating analysis sheets...")

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Sites', index=False)

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
        perform_random_lysine_analysis(df).to_excel(writer, sheet_name='Random_Lysine', index=False)

    print(f"Analysis complete. Output: {output_file}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
