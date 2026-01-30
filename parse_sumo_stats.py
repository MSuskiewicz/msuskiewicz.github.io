#!/usr/bin/env python3
"""
SUMO Site Analyzer (Stats Extended)
===================================
Analyzes SUMO sites using AlphaFold models.
- Categorizes sites: Consensus vs Acidic Flank vs Neither.
- Calculates Mean and Std Dev for SUMO scores per category.
- Output: Excel file with 'Sites' and statistical 'Analysis' sheets.

Usage:
    python parse_sumo_stats.py <input_excel_file> [output_file]
"""

import sys
import os
import re
import time
import math
import requests
import pandas as pd
import numpy as np
from scipy import stats
import logging

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

# =============================================================================

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

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
        r = requests.get(url, timeout=timeout)
        return r if r.status_code == 200 else None
    except Exception:
        return None

def calc_distance(c1: tuple, c2: tuple) -> float:
    return math.sqrt(sum((a - b)**2 for a, b in zip(c1, c2)))

# =============================================================================
# SEQUENCE & MAPPING LOGIC
# =============================================================================

def get_unisave_sequence(uid: str, cache: dict, min_length: int = 0) -> str | None:
    """
    Fetch sequence from UniSave (historical UniProt versions).
    Tries multiple versions including higher ones (up to 50) to find sequences
    that are long enough for the required position.

    Args:
        uid: UniProt accession
        cache: Cache dictionary
        min_length: Minimum required sequence length (0 = any length)
    """
    # If min_length specified, use a different cache key
    key = f"unisave_{uid}" if min_length == 0 else f"unisave_{uid}_min{min_length}"
    if key in cache:
        return cache[key]

    # Try a range of versions - early ones first, then jump to higher ones
    # Higher versions often have more complete/updated sequences
    versions_to_try = [1, 2, 3, 5, 10, 20, 30, 40, 50]

    best_seq = None
    for version in versions_to_try:
        r = safe_get(f"https://rest.uniprot.org/unisave/{uid}?format=fasta&versions={version}")
        if r and r.text.strip() and not r.text.startswith('<!'):
            lines = r.text.strip().split('\n')
            seq = ''.join(l for l in lines if not l.startswith('>'))
            if seq and len(seq) > 10:
                # If we need a minimum length, keep trying until we find one
                if min_length > 0:
                    if len(seq) >= min_length:
                        cache[key] = seq
                        return seq
                    # Keep track of the longest sequence found so far
                    if best_seq is None or len(seq) > len(best_seq):
                        best_seq = seq
                else:
                    # No minimum length requirement, return first valid sequence
                    cache[key] = seq
                    return seq

    # If we had a min_length requirement but couldn't meet it, return best we found
    if best_seq:
        cache[key] = best_seq
        return best_seq

    cache[key] = None
    return None

def get_uniprot_sequence(uid: str, isoform: int | None, cache: dict) -> str | None:
    key = f"upseq_{uid}_{isoform}"
    if key in cache: return cache[key]
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
    if key in cache: return cache[key]
    r = safe_get(f"https://rest.uniprot.org/uniprotkb/{uid}.json")
    if r:
        try:
            data = r.json()
            genes = [g['geneName']['value'] for g in data.get('genes', []) if 'geneName' in g]
            cache[key] = {'genes': genes, 'sequence': data.get('sequence', {}).get('value', ''), 'primary': data.get('primaryAccession', uid)}
        except Exception: cache[key] = None
    else: cache[key] = None
    return cache[key]

def get_mapped_accession(uid: str, cache: dict) -> str | None:
    key = f"mapped_{uid}"
    if key in cache: return cache[key]
    info = get_uniprot_info(uid, cache)
    if info and info['primary'] and info['primary'] != uid:
        cache[key] = info['primary']; return cache[key]
    r = safe_get(f"https://rest.uniprot.org/uniprotkb/search?query=(sec_acc:{uid})&fields=accession&size=1")
    if r:
        try:
            results = r.json().get('results', [])
            if results:
                acc = results[0].get('primaryAccession')
                if acc: cache[key] = acc; return acc
        except Exception: pass
    cache[key] = None
    return None

def sequence_similarity(s1: str, s2: str) -> float:
    if not s1 or not s2: return 0.0
    if len(s1) > len(s2): s1, s2 = s2, s1
    w = min(len(s1), 40)
    best = 0.0
    for i in range(0, max(1, len(s1) - w + 1), max(1, len(s1) // 4)):
        for j in range(len(s2) - w + 1):
            matches = sum(a == b for a, b in zip(s1[i:i+w], s2[j:j+w]))
            best = max(best, matches / w)
    return best

def map_position_by_residue(src_seq: str, tgt_seq: str, src_pos: int) -> int | None:
    if not src_seq or not tgt_seq or src_pos < 1 or src_pos > len(src_seq): return None
    target_aa = src_seq[src_pos - 1]
    for w in [20, 15, 10, 7, 5]:
        start, end = max(0, src_pos - 1 - w), min(len(src_seq), src_pos + w)
        query, pos_in_q = src_seq[start:end], src_pos - 1 - start
        best_pos, best_score = None, 0
        for i in range(max(0, len(tgt_seq) - len(query) + 1)):
            tgt_window = tgt_seq[i:i+len(query)]
            if len(tgt_window) < len(query): continue
            matches = sum(a == b for a, b in zip(query, tgt_window))
            tgt_idx = i + pos_in_q
            if tgt_idx < len(tgt_seq) and tgt_seq[tgt_idx] == target_aa:
                score = matches + 10
                if score > best_score: best_score, best_pos = score, tgt_idx + 1
        if best_pos and best_score >= len(query) * 0.4 + 10: return best_pos
    return None

# =============================================================================
# ALPHAFOLD HANDLING
# =============================================================================

def fetch_alphafold(uid: str, fragment: int = 1) -> dict | None:
    r = safe_get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}")
    if not r: return None
    try:
        data = r.json()
        if not data: return None
        cif_url = next((e['cifUrl'] for e in data if f"-F{fragment}-" in e.get('cifUrl', '')), None)
        if not cif_url and fragment == 1 and data: cif_url = data[0].get('cifUrl')
        if not cif_url: return None
        cif_r = safe_get(cif_url)
        return parse_cif(cif_r.text) if cif_r else None
    except Exception: return None

def parse_cif(content: str) -> dict:
    plddt, seq, coords, acidic = {}, {}, {}, []
    in_atom, headers, idx = False, [], {}
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('_atom_site.'):
            in_atom = True; name = line.split('.')[1].split()[0]
            idx[name] = len(headers); headers.append(name)
        elif in_atom and line.startswith(('ATOM', 'HETATM')):
            p = line.split()
            try:
                if p[idx['label_atom_id']] == 'CA':
                    rn = int(p[idx['label_seq_id']]); aa = AA_3TO1.get(p[idx['label_comp_id']], 'X')
                    plddt[rn] = float(p[idx['B_iso_or_equiv']]); seq[rn] = aa
                    coords[rn] = (float(p[idx['Cartn_x']]), float(p[idx['Cartn_y']]), float(p[idx['Cartn_z']]))
                    if aa in ACIDIC_RESIDUES: acidic.append(rn)
            except Exception: continue
        elif in_atom and line.startswith('#'): break
    seq_str = ''.join(seq.get(i, 'X') for i in range(1, max(seq.keys(), default=0) + 1))
    return {'plddt': plddt, 'sequence': seq, 'sequence_string': seq_str, 'coords': coords, 'acidic': acidic}

# =============================================================================
# RESOLUTION & ANALYSIS
# =============================================================================

def resolve_structure_and_position(uid: str, isoform: int | None, pos: int, cache: dict) -> tuple[dict | None, str, int]:
    cache_key = f"resolve_{uid}_{isoform}_{pos}"
    if cache_key in cache: return cache[cache_key]

    struct = fetch_alphafold(uid, 1)
    if struct and pos > max(struct['plddt'].keys(), default=0):
        for f in range(2, 21):
            next_struct = fetch_alphafold(uid, f)
            if not next_struct: break
            if pos in next_struct['plddt']: struct = next_struct; break

    if struct:
        final_pos = pos
        if isoform:
            iso_seq = get_uniprot_sequence(uid, isoform, cache)
            if iso_seq:
                mapped = map_position_by_residue(iso_seq, struct['sequence_string'], pos)
                if mapped: final_pos = mapped
                else: return (None, uid, pos)
        result = (struct, uid, final_pos)
        cache[cache_key] = result; return result

    # No direct AlphaFold model - try to find via mapping
    logger.debug(f"{uid}: No direct AlphaFold model, trying mapping...")

    # Try to get sequence with minimum length requirement matching the position
    original_seq = get_unisave_sequence(uid, cache, min_length=pos) or get_uniprot_sequence(uid, isoform, cache)
    if not original_seq or pos > len(original_seq):
        logger.warning(f"{uid}: No sequence found or position {pos} > length {len(original_seq) if original_seq else 0}")
        return (None, uid, pos)

    logger.debug(f"{uid}: Found sequence of length {len(original_seq)}")

    target_acc = get_mapped_accession(uid, cache)
    if target_acc:
        logger.debug(f"{uid}: Mapped to {target_acc}")
        struct = fetch_alphafold(target_acc, 1)
        if struct:
            mapped_pos = map_position_by_residue(original_seq, struct['sequence_string'], pos)
            if mapped_pos:
                logger.debug(f"{uid}: Position {pos} mapped to {mapped_pos} in {target_acc}")
                return (struct, target_acc, mapped_pos)
            else:
                logger.warning(f"{uid}: Failed to map position {pos} to {target_acc}")
        else:
            logger.warning(f"{uid}: Mapped acc {target_acc} has no AlphaFold model")
    else:
        logger.warning(f"{uid}: No mapped accession found")

    return (None, uid, pos)

def analyze_site(struct: dict, pos: int) -> dict:
    res = {k: None for k in ['plddt_site', 'plddt_11avg', 'aa_m2', 'aa_m1', 'aa_p1', 'aa_p2', 'dist_de_any', 'closest_de_any', 'dist_de_not_flank', 'closest_de_not_flank', 'forward_consensus', 'inverse_consensus', 'flexible', 'structured', 'de_in_space', 'de_in_flank', 'category']}
    if not struct or pos not in struct['plddt']: return res

    plddt, seq, coords, acidic = struct['plddt'], struct['sequence'], struct['coords'], struct['acidic']
    res['plddt_site'] = plddt.get(pos)
    window_vals = [plddt[i] for i in range(pos - 5, pos + 6) if i in plddt]
    res['plddt_11avg'] = sum(window_vals) / len(window_vals) if window_vals else None
    res['aa_m2'], res['aa_m1'], res['aa_p1'], res['aa_p2'] = seq.get(pos-2), seq.get(pos-1), seq.get(pos+1), seq.get(pos+2)

    if pos in coords and acidic:
        site_coord = coords[pos]
        min_any, cls_any = float('inf'), None
        min_nf, cls_nf = float('inf'), None
        for ac in acidic:
            d = calc_distance(site_coord, coords[ac])
            if d < min_any: min_any, cls_any = d, ac
            if abs(ac - pos) > 2 and d < min_nf: min_nf, cls_nf = d, ac
        if cls_any: res['dist_de_any'], res['closest_de_any'] = round(min_any, 2), f"{seq.get(cls_any, '?')}{cls_any}"
        if cls_nf: res['dist_de_not_flank'], res['closest_de_not_flank'] = round(min_nf, 2), f"{seq.get(cls_nf, '?')}{cls_nf}"

    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    res['inverse_consensus'] = 'Yes' if res['aa_p1'] in HYDROPHOBIC_RESIDUES and res['aa_m2'] in ACIDIC_RESIDUES else None

    if res['plddt_11avg'] is not None:
        if res['plddt_11avg'] < PLDDT_THRESHOLD: res['flexible'] = 'Yes'
        else: res['structured'] = 'Yes'

    if res['dist_de_not_flank'] is not None and res['dist_de_not_flank'] <= DISTANCE_THRESHOLD: res['de_in_space'] = 'Yes'
    if res['aa_m2'] in ACIDIC_RESIDUES or res['aa_p2'] in ACIDIC_RESIDUES: res['de_in_flank'] = 'Yes'

    is_consensus = (res['forward_consensus'] == 'Yes' or res['inverse_consensus'] == 'Yes')
    is_acidic_flank = (res['de_in_flank'] == 'Yes')

    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"

    if is_consensus: suffix = "_consensus"
    elif is_acidic_flank: suffix = "_acidic_flank"
    else: suffix = "_neither"

    res['category'] = prefix + suffix
    return res

# =============================================================================
# STATISTICS & MAIN
# =============================================================================

def perform_analysis(df: pd.DataFrame, score_category: str) -> pd.DataFrame:
    score_col = 'Score (SUMO site)'
    data_list = []

    cats = [
        'Flexible_consensus',
        'Flexible_acidic_flank',
        'Flexible_neither',
        'Structured_consensus',
        'Structured_acidic_flank',
        'Structured_neither'
    ]

    # 1. Category Statistics (Count, Mean Score, Std Score)
    for cat in cats:
        subset = df[df['Category'] == cat]
        count = len(subset)
        if count > 0 and score_col in df.columns:
            mean_score = subset[score_col].mean()
            std_score = subset[score_col].std()
        else:
            mean_score = 0.0
            std_score = 0.0

        data_list.append({
            'Metric': cat,
            'Value': count,
            'Mean_Score': round(mean_score, 2) if not pd.isna(mean_score) else None,
            'Std_Score': round(std_score, 2) if not pd.isna(std_score) else None
        })

    # Spacer
    data_list.append({'Metric': '', 'Value': '', 'Mean_Score': '', 'Std_Score': ''})

    # 2. Consensus Rates
    total_flex = len(df[df['Flexible'] == 'Yes'])
    total_struct = len(df[df['Structured'] == 'Yes'])

    flex_cons = len(df[df['Category'] == 'Flexible_consensus'])
    struct_cons = len(df[df['Category'] == 'Structured_consensus'])

    f_rate = flex_cons / total_flex if total_flex > 0 else 0
    s_rate = struct_cons / total_struct if total_struct > 0 else 0

    data_list.append({
        'Metric': 'Flexible consensus rate',
        'Value': f"{f_rate:.3f} ({flex_cons}/{total_flex})",
        'Mean_Score': None, 'Std_Score': None
    })

    data_list.append({
        'Metric': 'Structured consensus rate',
        'Value': f"{s_rate:.3f} ({struct_cons}/{total_struct})",
        'Mean_Score': None, 'Std_Score': None
    })

    # 3. Chi-Square Test (2x3: Flexible vs Structured across the 3 sub-types)
    f_counts = [len(df[df['Category'] == c]) for c in cats[:3]]
    s_counts = [len(df[df['Category'] == c]) for c in cats[3:]]

    if sum(f_counts) > 0 and sum(s_counts) > 0:
        try:
            chi2, pval, _, _ = stats.chi2_contingency([f_counts, s_counts])
            data_list.append({'Metric': 'Chi-square (2x3)', 'Value': f"{chi2:.2f}", 'Mean_Score': None, 'Std_Score': None})
            data_list.append({'Metric': 'p-value', 'Value': f"{pval:.2e}", 'Mean_Score': None, 'Std_Score': None})
            data_list.append({'Metric': 'Conclusion', 'Value': 'Significant' if pval < 0.05 else 'Not Significant', 'Mean_Score': None, 'Std_Score': None})
        except Exception as e:
            data_list.append({'Metric': 'Chi-square error', 'Value': str(e), 'Mean_Score': None, 'Std_Score': None})

    # Return DataFrame with specific column order
    return pd.DataFrame(data_list, columns=['Metric', 'Value', 'Mean_Score', 'Std_Score'])

def process_file(input_file: str, output_file: str = None):
    print(f"\n{'='*60}\nSUMO Site Analyzer (Stats Extended)\n{'='*60}\n")
    df = pd.read_excel(input_file, sheet_name=0, header=1)
    protein_col = next((c for c in df.columns if str(c).lower() == 'protein'), 'Protein')
    position_col = next((c for c in df.columns if str(c).lower() == 'position'), 'Position')

    cache, results = {}, []
    for idx, row in df.iterrows():
        if (idx + 1) % 50 == 0: logger.info(f"Processing {idx + 1}/{len(df)}")
        uid, iso = parse_uniprot_id(str(row[protein_col]).strip())
        try: pos = int(row[position_col])
        except: results.append({'used_id': uid, 'mapped_pos': None}); continue

        struct, used_id, final_pos = resolve_structure_and_position(uid, iso, pos, cache)
        analysis = analyze_site(struct, final_pos)
        analysis.update({'used_id': used_id, 'mapped_pos': final_pos if final_pos != pos else None})
        results.append(analysis); time.sleep(0.01)

    res_df = pd.DataFrame(results)
    cols_map = {
        'used_id': 'UniProt_ID_used', 'mapped_pos': 'Mapped_position',
        'plddt_site': 'pLDDT_site', 'plddt_11avg': 'pLDDT_11residue_avg',
        'aa_m2': 'AA_minus2', 'aa_m1': 'AA_minus1', 'aa_p1': 'AA_plus1', 'aa_p2': 'AA_plus2',
        'dist_de_any': 'Dist_DE_any_A', 'closest_de_any': 'Closest_DE_any',
        'dist_de_not_flank': 'Dist_DE_not_flank_A', 'closest_de_not_flank': 'Closest_DE_not_flank',
        'forward_consensus': 'Forward_consensus', 'inverse_consensus': 'Inverse_consensus',
        'flexible': 'Flexible', 'structured': 'Structured',
        'de_in_space': 'Asp/Glu_in_space', 'de_in_flank': 'Asp/Glu_in_-2/+2', 'category': 'Category'
    }
    res_df = res_df.rename(columns=cols_map)
    df = pd.concat([df, res_df], axis=1)

    # Score stratification & Type Conversion
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

    if not output_file: output_file = os.path.splitext(input_file)[0] + '_analyzed_stats.xlsx'
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Sites', index=False)
        for name, sub in cats.items():
            if len(sub) > 0: perform_analysis(sub, name).to_excel(writer, sheet_name=f'Analysis_{name}', index=False)
    print(f"Analysis complete. Output: {output_file}")

def main():
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)

if __name__ == "__main__": main()
