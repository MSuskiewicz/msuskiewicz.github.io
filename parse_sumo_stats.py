#!/usr/bin/env python3
"""
SUMO Site Analyzer (Stats Extended)
===================================
Analyzes SUMO sites using AlphaFold models.
- Categorizes sites: Consensus, Flank_acidic, Space_acidic, Neither.
- Measures distances using NZ (Lys) to CG (Asp) / CD (Glu).
- Calculates Mean and Std Dev for SUMO scores per category.
- Output: Excel file with 'Sites', statistical 'Analysis' sheets, and trend analysis.

Usage:
    python parse_sumo_stats.py <input_excel_file> [output_file]
"""

import sys
import os
import re
import math
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
# ALPHAFOLD HANDLING - Now extracts NZ (Lys) and CG/CD (acidic) coordinates
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
    """
    Parse mmCIF file to extract:
    - pLDDT (B-factor) per residue
    - Sequence
    - CA coordinates for all residues
    - NZ coordinates for lysine residues
    - CG coordinates for Asp, CD coordinates for Glu (acidic sidechain atoms)
    """
    plddt, seq, coords_ca, coords_nz, acidic_coords = {}, {}, {}, {}, {}
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

                # CA atom: store pLDDT, sequence, and CA coordinates
                if atom_name == 'CA':
                    plddt[rn] = float(p[idx['B_iso_or_equiv']])
                    seq[rn] = aa
                    coords_ca[rn] = coord

                # NZ atom of lysine
                elif atom_name == 'NZ' and aa3 == 'LYS':
                    coords_nz[rn] = coord

                # CG atom of Asp (sidechain carboxyl)
                elif atom_name == 'CG' and aa3 == 'ASP':
                    acidic_coords[rn] = coord

                # CD atom of Glu (sidechain carboxyl)
                elif atom_name == 'CD' and aa3 == 'GLU':
                    acidic_coords[rn] = coord

            except Exception:
                continue
        elif in_atom and line.startswith('#'):
            break

    seq_str = ''.join(seq.get(i, 'X') for i in range(1, max(seq.keys(), default=0) + 1))

    # List of acidic residue positions
    acidic_positions = sorted(acidic_coords.keys())

    return {
        'plddt': plddt,
        'sequence': seq,
        'sequence_string': seq_str,
        'coords_ca': coords_ca,
        'coords_nz': coords_nz,
        'acidic_coords': acidic_coords,
        'acidic_positions': acidic_positions
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

    # No direct AlphaFold model - try to find via mapping
    original_seq = get_unisave_sequence(uid, cache, min_length=pos) or get_uniprot_sequence(uid, isoform, cache)
    if not original_seq or pos > len(original_seq):
        return (None, uid, pos)

    target_acc = get_mapped_accession(uid, cache)
    if target_acc:
        struct = fetch_alphafold(target_acc, 1)
        if struct:
            mapped_pos = map_position_by_residue(original_seq, struct['sequence_string'], pos)
            if mapped_pos:
                return (struct, target_acc, mapped_pos)

    return (None, uid, pos)


def analyze_site(struct: dict, pos: int) -> dict:
    """
    Analyze a SUMO site for:
    - pLDDT values (site and 11-residue average)
    - Flanking amino acids
    - Distance to acidic residues (using NZ of Lys to CG/CD of acidic)
    - Consensus motifs (forward: ψKxE, inverse: ExKψ)
    - Flexibility classification
    - Category assignment: consensus, flank_acidic, space_acidic, neither
    """
    res = {
        'plddt_site': None, 'plddt_11avg': None,
        'aa_m2': None, 'aa_m1': None, 'aa_p1': None, 'aa_p2': None,
        'dist_acidic_any': None, 'closest_acidic_any': None,
        'dist_acidic_space': None, 'closest_acidic_space': None,
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

        for ac in acidic_positions:
            if ac not in acidic_coords:
                continue
            d = calc_distance(lys_nz_coord, acidic_coords[ac])

            # Track closest acidic (any position)
            if d < min_any:
                min_any, cls_any = d, ac

            # Track closest acidic NOT at position -2 or +2 (for "in space" calculation)
            # Note: positions -1, +1, and all other positions are allowed
            rel_pos = ac - pos
            if rel_pos != -2 and rel_pos != 2:
                if d < min_space:
                    min_space, cls_space = d, ac

        if cls_any:
            res['dist_acidic_any'] = round(min_any, 2)
            res['closest_acidic_any'] = f"{seq.get(cls_any, '?')}{cls_any}"
        if cls_space:
            res['dist_acidic_space'] = round(min_space, 2)
            res['closest_acidic_space'] = f"{seq.get(cls_space, '?')}{cls_space}"

    # Consensus motifs
    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    res['inverse_consensus'] = 'Yes' if res['aa_p1'] in HYDROPHOBIC_RESIDUES and res['aa_m2'] in ACIDIC_RESIDUES else None

    # Flexibility classification based on 11-residue average pLDDT
    if res['plddt_11avg'] is not None:
        if res['plddt_11avg'] < PLDDT_THRESHOLD:
            res['flexible'] = 'Yes'
        else:
            res['structured'] = 'Yes'

    # Acidic in space: within distance threshold, excluding positions -2 and +2
    if res['dist_acidic_space'] is not None and res['dist_acidic_space'] <= DISTANCE_THRESHOLD:
        res['acidic_in_space'] = 'Yes'

    # Acidic in flank: specifically at position -2 or +2
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
# STATISTICS & ANALYSIS
# =============================================================================

def perform_analysis(df: pd.DataFrame, score_category: str, summary_stats: dict = None) -> pd.DataFrame:
    """Perform statistical analysis on the sites data."""
    score_col = 'Score (SUMO site)'
    data_list = []

    # Data completeness summary
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

    # Category counts
    data_list.append({'Metric': '=== CATEGORY COUNTS ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Flexible categories
    add_stat_row('Flexible_consensus', df[df['Category'] == 'Flexible_consensus'])
    add_stat_row('Flexible_flank_acidic', df[df['Category'] == 'Flexible_flank_acidic'])
    add_stat_row('Flexible_space_acidic', df[df['Category'] == 'Flexible_space_acidic'])
    add_stat_row('Flexible_neither', df[df['Category'] == 'Flexible_neither'])

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Structured categories
    add_stat_row('Structured_consensus', df[df['Category'] == 'Structured_consensus'])
    add_stat_row('Structured_flank_acidic', df[df['Category'] == 'Structured_flank_acidic'])
    add_stat_row('Structured_space_acidic', df[df['Category'] == 'Structured_space_acidic'])
    add_stat_row('Structured_neither', df[df['Category'] == 'Structured_neither'])

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Rates - cumulative categories
    data_list.append({'Metric': '=== RATES (cumulative) ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    total_flex = len(df[df['Flexible'] == 'Yes'])
    total_struct = len(df[df['Structured'] == 'Yes'])

    # Count per category
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

    # Consensus rate
    data_list.append({'Metric': 'Flexible consensus rate', 'Value': rate_str(flex_cons, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured consensus rate', 'Value': rate_str(struct_cons, total_struct), 'Mean_Score': None, 'Std_Score': None})

    # Consensus + flank_acidic rate
    flex_cons_flank = flex_cons + flex_flank
    struct_cons_flank = struct_cons + struct_flank
    data_list.append({'Metric': 'Flexible consensus+flank rate', 'Value': rate_str(flex_cons_flank, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured consensus+flank rate', 'Value': rate_str(struct_cons_flank, total_struct), 'Mean_Score': None, 'Std_Score': None})

    # All acidic rate (consensus + flank + space)
    flex_all_acidic = flex_cons + flex_flank + flex_space
    struct_all_acidic = struct_cons + struct_flank + struct_space
    data_list.append({'Metric': 'Flexible all_acidic rate', 'Value': rate_str(flex_all_acidic, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured all_acidic rate', 'Value': rate_str(struct_all_acidic, total_struct), 'Mean_Score': None, 'Std_Score': None})

    # None rate
    data_list.append({'Metric': 'Flexible none rate', 'Value': rate_str(flex_neither, total_flex), 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured none rate', 'Value': rate_str(struct_neither, total_struct), 'Mean_Score': None, 'Std_Score': None})

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Chi-Square Test: Flexible vs Structured across 4 categories
    data_list.append({'Metric': '=== STATISTICAL TESTS ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    data_list.append({'Metric': 'Chi-square: Flex vs Struct distribution', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    data_list.append({'Metric': '(across consensus/flank/space/neither)', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    f_counts = [flex_cons, flex_flank, flex_space, flex_neither]
    s_counts = [struct_cons, struct_flank, struct_space, struct_neither]

    if sum(f_counts) > 0 and sum(s_counts) > 0:
        try:
            chi2, pval, _, _ = stats.chi2_contingency([f_counts, s_counts])
            data_list.append({'Metric': 'Chi-square statistic', 'Value': f"{chi2:.2f}", 'Mean_Score': None, 'Std_Score': None})
            data_list.append({'Metric': 'p-value', 'Value': f"{pval:.2e}", 'Mean_Score': None, 'Std_Score': None})
            data_list.append({'Metric': 'Conclusion', 'Value': 'Significant difference' if pval < 0.05 else 'No significant difference', 'Mean_Score': None, 'Std_Score': None})
        except Exception as e:
            data_list.append({'Metric': 'Chi-square error', 'Value': str(e), 'Mean_Score': None, 'Std_Score': None})

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Inverse correlation test: Acidic_in_-2/+2 vs Acidic_in_space
    data_list.append({'Metric': '=== INVERSE CORRELATION TEST ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    data_list.append({'Metric': 'Testing: Acidic_in_flank vs Acidic_in_space', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Build 2x2 contingency table
    flank_yes_space_yes = len(df[(df['Acidic_in_-2/+2'] == 'Yes') & (df['Acidic_in_space'] == 'Yes')])
    flank_yes_space_no = len(df[(df['Acidic_in_-2/+2'] == 'Yes') & (df['Acidic_in_space'] != 'Yes')])
    flank_no_space_yes = len(df[(df['Acidic_in_-2/+2'] != 'Yes') & (df['Acidic_in_space'] == 'Yes')])
    flank_no_space_no = len(df[(df['Acidic_in_-2/+2'] != 'Yes') & (df['Acidic_in_space'] != 'Yes')])

    contingency_flank_space = [[flank_yes_space_yes, flank_yes_space_no], [flank_no_space_yes, flank_no_space_no]]

    data_list.append({'Metric': 'Flank=Yes & Space=Yes', 'Value': flank_yes_space_yes, 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Flank=Yes & Space=No', 'Value': flank_yes_space_no, 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Flank=No & Space=Yes', 'Value': flank_no_space_yes, 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Flank=No & Space=No', 'Value': flank_no_space_no, 'Mean_Score': None, 'Std_Score': None})

    try:
        odds_ratio, fisher_p = stats.fisher_exact(contingency_flank_space)
        data_list.append({'Metric': 'Fisher exact odds ratio', 'Value': f"{odds_ratio:.3f}", 'Mean_Score': None, 'Std_Score': None})
        data_list.append({'Metric': 'Fisher exact p-value', 'Value': f"{fisher_p:.2e}", 'Mean_Score': None, 'Std_Score': None})
        if fisher_p < 0.05:
            if odds_ratio < 1:
                conclusion = 'Significant INVERSE correlation (mutually exclusive)'
            else:
                conclusion = 'Significant POSITIVE correlation (co-occurring)'
        else:
            conclusion = 'No significant correlation'
        data_list.append({'Metric': 'Conclusion', 'Value': conclusion, 'Mean_Score': None, 'Std_Score': None})
    except Exception as e:
        data_list.append({'Metric': 'Fisher test error', 'Value': str(e), 'Mean_Score': None, 'Std_Score': None})

    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Test: Do Structured sites without consensus/flank have enriched Acidic_in_space?
    data_list.append({'Metric': '=== SPACE ACIDIC ENRICHMENT TEST ===', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})
    data_list.append({'Metric': 'Testing: Structured_neither has more Acidic_in_space?', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # Structured sites with consensus or flank_acidic
    struct_with_cons_flank = df[(df['Structured'] == 'Yes') & (df['Category'].isin(['Structured_consensus', 'Structured_flank_acidic']))]
    struct_with_cons_flank_space_yes = len(struct_with_cons_flank[struct_with_cons_flank['Acidic_in_space'] == 'Yes'])
    struct_with_cons_flank_total = len(struct_with_cons_flank)

    # Structured sites with neither (no consensus, no flank_acidic, but may have space_acidic)
    struct_neither_df = df[(df['Structured'] == 'Yes') & (~df['Category'].isin(['Structured_consensus', 'Structured_flank_acidic']))]
    struct_neither_space_yes = len(struct_neither_df[struct_neither_df['Acidic_in_space'] == 'Yes'])
    struct_neither_total = len(struct_neither_df)

    rate_cons_flank = struct_with_cons_flank_space_yes / struct_with_cons_flank_total if struct_with_cons_flank_total > 0 else 0
    rate_neither = struct_neither_space_yes / struct_neither_total if struct_neither_total > 0 else 0

    data_list.append({'Metric': 'Structured (cons/flank) Acidic_in_space rate', 'Value': f"{rate_cons_flank:.3f} ({struct_with_cons_flank_space_yes}/{struct_with_cons_flank_total})", 'Mean_Score': None, 'Std_Score': None})
    data_list.append({'Metric': 'Structured (other) Acidic_in_space rate', 'Value': f"{rate_neither:.3f} ({struct_neither_space_yes}/{struct_neither_total})", 'Mean_Score': None, 'Std_Score': None})

    # 2x2 Fisher test
    contingency_enrichment = [
        [struct_with_cons_flank_space_yes, struct_with_cons_flank_total - struct_with_cons_flank_space_yes],
        [struct_neither_space_yes, struct_neither_total - struct_neither_space_yes]
    ]

    try:
        if struct_with_cons_flank_total > 0 and struct_neither_total > 0:
            odds_ratio, fisher_p = stats.fisher_exact(contingency_enrichment)
            data_list.append({'Metric': 'Fisher exact odds ratio', 'Value': f"{odds_ratio:.3f}", 'Mean_Score': None, 'Std_Score': None})
            data_list.append({'Metric': 'Fisher exact p-value', 'Value': f"{fisher_p:.2e}", 'Mean_Score': None, 'Std_Score': None})
            if fisher_p < 0.05 and rate_neither > rate_cons_flank:
                conclusion = 'Significant enrichment in non-consensus structured sites'
            elif fisher_p < 0.05:
                conclusion = 'Significant but lower in non-consensus structured sites'
            else:
                conclusion = 'No significant difference'
            data_list.append({'Metric': 'Conclusion', 'Value': conclusion, 'Mean_Score': None, 'Std_Score': None})
    except Exception as e:
        data_list.append({'Metric': 'Enrichment test error', 'Value': str(e), 'Mean_Score': None, 'Std_Score': None})

    return pd.DataFrame(data_list, columns=['Metric', 'Value', 'Mean_Score', 'Std_Score'])


def perform_trend_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """
    Global trend analysis:
    1. Does higher score correlate with larger difference between Flexible
       and Structured in acidic rates?
    2. Does score correlate with presence of any acidic nearby (flank or space)?

    Uses Spearman correlation for analysis.
    """
    score_col = 'Score (SUMO site)'
    data_list = []

    data_list.append({'Metric': '=== GLOBAL SCORE-ACIDIC CORRELATION ===', 'Value': ''})
    data_list.append({'Metric': 'Testing: Does score correlate with acidic presence?', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    if score_col not in df.columns:
        data_list.append({'Metric': 'Error', 'Value': 'Score column not found'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    # Filter to rows with valid pLDDT and score
    valid_df = df[df['pLDDT_site'].notna() & df[score_col].notna()].copy()

    # First: Global correlation between score and acidic presence
    if len(valid_df) > 10:
        # Create binary columns for acidic presence
        valid_df['has_any_acidic'] = ((valid_df['Acidic_in_-2/+2'] == 'Yes') |
                                       (valid_df['Acidic_in_space'] == 'Yes')).astype(int)
        valid_df['has_flank_acidic'] = (valid_df['Acidic_in_-2/+2'] == 'Yes').astype(int)
        valid_df['has_space_acidic'] = (valid_df['Acidic_in_space'] == 'Yes').astype(int)
        valid_df['is_consensus'] = ((valid_df['Forward_consensus'] == 'Yes') |
                                     (valid_df['Inverse_consensus'] == 'Yes')).astype(int)

        # Point-biserial correlation (using Spearman for binary variable)
        rho_any, p_any = stats.spearmanr(valid_df[score_col], valid_df['has_any_acidic'])
        rho_flank, p_flank = stats.spearmanr(valid_df[score_col], valid_df['has_flank_acidic'])
        rho_space, p_space = stats.spearmanr(valid_df[score_col], valid_df['has_space_acidic'])
        rho_cons, p_cons = stats.spearmanr(valid_df[score_col], valid_df['is_consensus'])

        data_list.append({'Metric': 'Score vs any acidic (flank or space)', 'Value': f"rho={rho_any:.3f}, p={p_any:.2e}"})
        data_list.append({'Metric': 'Score vs flank acidic (-2/+2)', 'Value': f"rho={rho_flank:.3f}, p={p_flank:.2e}"})
        data_list.append({'Metric': 'Score vs space acidic', 'Value': f"rho={rho_space:.3f}, p={p_space:.2e}"})
        data_list.append({'Metric': 'Score vs consensus', 'Value': f"rho={rho_cons:.3f}, p={p_cons:.2e}"})

        # Interpretation
        if p_any < 0.05:
            if rho_any > 0:
                data_list.append({'Metric': 'Conclusion', 'Value': 'Higher score = more likely to have nearby acidic'})
            else:
                data_list.append({'Metric': 'Conclusion', 'Value': 'Higher score = less likely to have nearby acidic'})
        else:
            data_list.append({'Metric': 'Conclusion', 'Value': 'No significant correlation between score and acidic presence'})

        data_list.append({'Metric': '', 'Value': ''})

        # Also test separately for Flexible vs Structured
        data_list.append({'Metric': '--- By flexibility ---', 'Value': ''})

        flex_df = valid_df[valid_df['Flexible'] == 'Yes']
        struct_df = valid_df[valid_df['Structured'] == 'Yes']

        if len(flex_df) > 10:
            rho_f, p_f = stats.spearmanr(flex_df[score_col], flex_df['has_any_acidic'])
            data_list.append({'Metric': 'Flexible: Score vs any acidic', 'Value': f"rho={rho_f:.3f}, p={p_f:.2e}"})

        if len(struct_df) > 10:
            rho_s, p_s = stats.spearmanr(struct_df[score_col], struct_df['has_any_acidic'])
            data_list.append({'Metric': 'Structured: Score vs any acidic', 'Value': f"rho={rho_s:.3f}, p={p_s:.2e}"})

        data_list.append({'Metric': '', 'Value': ''})

    data_list.append({'Metric': '=== FLEX-STRUCT DIFFERENCE TREND ===', 'Value': ''})
    data_list.append({'Metric': 'Testing: Does higher score predict larger Flex-Struct difference?', 'Value': ''})
    data_list.append({'Metric': '', 'Value': ''})

    if len(valid_df) < 50:
        data_list.append({'Metric': 'Warning', 'Value': f'Too few valid sites ({len(valid_df)}) for trend analysis'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    # Sort by score and create bins
    valid_df = valid_df.sort_values(score_col)
    n_bins = 10
    valid_df['score_bin'] = pd.qcut(valid_df[score_col], q=n_bins, labels=False, duplicates='drop')

    bin_results = []

    for bin_idx in sorted(valid_df['score_bin'].unique()):
        bin_df = valid_df[valid_df['score_bin'] == bin_idx]
        mean_score = bin_df[score_col].mean()

        # Calculate rates for this bin
        flex_sites = bin_df[bin_df['Flexible'] == 'Yes']
        struct_sites = bin_df[bin_df['Structured'] == 'Yes']

        total_flex = len(flex_sites)
        total_struct = len(struct_sites)

        if total_flex == 0 or total_struct == 0:
            continue

        # Consensus + flank rate
        flex_cons_flank = len(flex_sites[flex_sites['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic'])])
        struct_cons_flank = len(struct_sites[struct_sites['Category'].isin(['Structured_consensus', 'Structured_flank_acidic'])])

        flex_rate_cf = flex_cons_flank / total_flex
        struct_rate_cf = struct_cons_flank / total_struct
        diff_cf = flex_rate_cf - struct_rate_cf

        # All acidic rate
        flex_all_acidic = len(flex_sites[flex_sites['Category'].isin(['Flexible_consensus', 'Flexible_flank_acidic', 'Flexible_space_acidic'])])
        struct_all_acidic = len(struct_sites[struct_sites['Category'].isin(['Structured_consensus', 'Structured_flank_acidic', 'Structured_space_acidic'])])

        flex_rate_all = flex_all_acidic / total_flex
        struct_rate_all = struct_all_acidic / total_struct
        diff_all = flex_rate_all - struct_rate_all

        bin_results.append({
            'bin': bin_idx,
            'mean_score': mean_score,
            'n_sites': len(bin_df),
            'n_flex': total_flex,
            'n_struct': total_struct,
            'flex_rate_cf': flex_rate_cf,
            'struct_rate_cf': struct_rate_cf,
            'diff_cf': diff_cf,
            'flex_rate_all': flex_rate_all,
            'struct_rate_all': struct_rate_all,
            'diff_all': diff_all
        })

    if len(bin_results) < 3:
        data_list.append({'Metric': 'Warning', 'Value': 'Too few bins with both Flex and Struct sites'})
        return pd.DataFrame(data_list, columns=['Metric', 'Value'])

    bin_df = pd.DataFrame(bin_results)

    # Spearman correlation: score vs difference
    data_list.append({'Metric': '--- Consensus + Flank Analysis ---', 'Value': ''})

    rho_cf, p_cf = stats.spearmanr(bin_df['mean_score'], bin_df['diff_cf'])
    data_list.append({'Metric': 'Spearman rho (score vs Flex-Struct diff)', 'Value': f"{rho_cf:.3f}"})
    data_list.append({'Metric': 'p-value', 'Value': f"{p_cf:.2e}"})

    if p_cf < 0.05:
        if rho_cf > 0:
            conclusion = 'Higher score = larger Flex-Struct difference (Flex more acidic)'
        else:
            conclusion = 'Higher score = smaller Flex-Struct difference'
    else:
        conclusion = 'No significant trend'
    data_list.append({'Metric': 'Conclusion', 'Value': conclusion})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '--- All Acidic Analysis ---', 'Value': ''})

    rho_all, p_all = stats.spearmanr(bin_df['mean_score'], bin_df['diff_all'])
    data_list.append({'Metric': 'Spearman rho (score vs Flex-Struct diff)', 'Value': f"{rho_all:.3f}"})
    data_list.append({'Metric': 'p-value', 'Value': f"{p_all:.2e}"})

    if p_all < 0.05:
        if rho_all > 0:
            conclusion = 'Higher score = larger Flex-Struct difference (Flex more acidic)'
        else:
            conclusion = 'Higher score = smaller Flex-Struct difference'
    else:
        conclusion = 'No significant trend'
    data_list.append({'Metric': 'Conclusion', 'Value': conclusion})

    data_list.append({'Metric': '', 'Value': ''})
    data_list.append({'Metric': '=== BIN DETAILS ===', 'Value': ''})

    for r in bin_results:
        data_list.append({
            'Metric': f"Bin {r['bin']} (mean score={r['mean_score']:.0f}, n={r['n_sites']})",
            'Value': f"Flex CF={r['flex_rate_cf']:.2f}, Struct CF={r['struct_rate_cf']:.2f}, Diff={r['diff_cf']:.2f}"
        })

    return pd.DataFrame(data_list, columns=['Metric', 'Value'])


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_file(input_file: str, output_file: str = None):
    print(f"\n{'=' * 60}\nSUMO Site Analyzer (Stats Extended)\n{'=' * 60}\n")
    df = pd.read_excel(input_file, sheet_name=0, header=1)

    # Only keep columns A to AA (first 27 columns) from original input
    df = df.iloc[:, :27]

    protein_col = next((c for c in df.columns if str(c).lower() == 'protein'), 'Protein')
    position_col = next((c for c in df.columns if str(c).lower() == 'position'), 'Position')

    # Shared cache
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

    # Process sites in parallel
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
        'forward_consensus': 'Forward_consensus',
        'inverse_consensus': 'Inverse_consensus',
        'flexible': 'Flexible',
        'structured': 'Structured',
        'acidic_in_space': 'Acidic_in_space',
        'acidic_in_flank': 'Acidic_in_-2/+2',
        'category': 'Category'
    }
    res_df = res_df.rename(columns=cols_map)

    # Reorder columns: UniProt_ID_used and Mapped_position first, then the rest
    ordered_cols = [
        'UniProt_ID_used', 'Mapped_position', 'pLDDT_site', 'pLDDT_11residue_avg',
        'AA_minus2', 'AA_minus1', 'AA_plus1', 'AA_plus2',
        'Dist_acidic_any_A', 'Closest_acidic_any', 'Dist_acidic_space_A', 'Closest_acidic_space',
        'Forward_consensus', 'Inverse_consensus', 'Flexible', 'Structured',
        'Acidic_in_space', 'Acidic_in_-2/+2', 'Category'
    ]
    res_df = res_df[[c for c in ordered_cols if c in res_df.columns]]

    df = pd.concat([df, res_df], axis=1)

    # Calculate success/failure summary
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

    # Score stratification
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

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Sites', index=False)

        for name, sub in cats.items():
            if len(sub) > 0:
                cat_total = len(sub)
                cat_success = sub['pLDDT_site'].notna().sum()
                cat_failed = cat_total - cat_success
                cat_summary = {'total': cat_total, 'success': cat_success, 'failed': cat_failed}
                perform_analysis(sub, name, cat_summary).to_excel(writer, sheet_name=f'Analysis_{name}', index=False)

        # Global trend analysis sheet
        trend_df = perform_trend_analysis(df)
        trend_df.to_excel(writer, sheet_name='Trend_Analysis', index=False)

    print(f"Analysis complete. Output: {output_file}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
