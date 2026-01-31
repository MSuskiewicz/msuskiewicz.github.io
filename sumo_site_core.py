#!/usr/bin/env python3
"""
SUMO Site Analyzer - Core Output
================================
Parses Excel files with SUMO modification sites, fetches AlphaFold structures,
and generates core analysis data.

Usage:
    python sumo_site_core.py <input_excel_file> [output_file]
"""

import sys
import os
import re
import math
import requests
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

PLDDT_THRESHOLD = 65.0          # pLDDT below this = Flexible, above = Structured
PLDDT_WINDOW_LENGTH = 11        # Window size for pLDDT averaging (centered on site)
DISTANCE_THRESHOLD = 7.0        # Angstroms for acidic proximity

# Score thresholds for stratification
SCORE_VERY_HIGH_THRESHOLD = 600
SCORE_HIGH_THRESHOLD = 400
SCORE_MEDIUM_THRESHOLD = 200

# Sliding window for fraction calculations
SLIDING_WINDOW_LENGTH = 100  # Number of sites for sliding window fractions

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


def parse_uniprot_id(protein_id: str) -> tuple[str, int | None]:
    """Extract UniProt accession and isoform number from protein ID."""
    match = re.search(r'([A-Z][0-9][A-Z0-9]{3,9}[0-9])(?:-(\d+))?', protein_id)
    if match:
        return match.group(1), int(match.group(2)) if match.group(2) else None
    return protein_id, None


def safe_get(url: str, timeout: int = 30) -> requests.Response | None:
    """Safely perform HTTP GET request with error handling."""
    try:
        r = get_session().get(url, timeout=timeout)
        return r if r.status_code == 200 else None
    except Exception:
        return None


def calc_distance(c1: tuple, c2: tuple) -> float:
    """Calculate Euclidean distance between two 3D coordinates."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


# =============================================================================
# SEQUENCE & MAPPING LOGIC
# =============================================================================

def get_unisave_sequence(uid: str, cache: dict, min_length: int = 0) -> str | None:
    """Retrieve sequence from UniSave (historical UniProt versions)."""
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
    """Retrieve current UniProt sequence."""
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
    """Retrieve UniProt entry information."""
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
    """Get mapped/current accession for obsolete UniProt IDs."""
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
    """Map position from source sequence to target sequence using local alignment."""
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
    """Fetch AlphaFold structure for a UniProt accession."""
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
    """Parse mmCIF file to extract coordinates and pLDDT values.

    Extracts:
    - CA coordinates for all residues
    - NZ coordinates for lysines (for distance calculations)
    - CG coordinates for Asp, CD coordinates for Glu (acidic residues)
    - pLDDT values from B-factor column
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

                if atom_name == 'CA':
                    plddt[rn] = float(p[idx['B_iso_or_equiv']])
                    seq[rn] = aa
                    coords_ca[rn] = coord
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
        'acidic_positions': acidic_positions
    }


# =============================================================================
# RESOLUTION & ANALYSIS
# =============================================================================

def resolve_structure_and_position(uid: str, isoform: int | None, pos: int, cache: dict) -> tuple[dict | None, str, int]:
    """Resolve AlphaFold structure and map position if needed."""
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
    """Analyze a SUMO site and compute all metrics.

    Returns dict with all output columns.
    """
    res = {
        'plddt_site': None,
        'plddt_window_avg': None,
        'aa_m2': None,
        'aa_m1': None,
        'aa_site': None,
        'aa_p1': None,
        'aa_p2': None,
        'dist_acidic_any': None,
        'closest_acidic_any': None,
        'dist_acidic_not_flank': None,
        'closest_acidic_not_flank': None,
        'n_acidic_within_threshold': 0,
        'flexible': None,
        'structured': None,
        'forward_consensus': None,
        'inverse_consensus': None,
        'acidic_in_flank': None,
        'acidic_in_space': None,
        'any_acidic_within_threshold': None,
        'category': None
    }

    if not struct or pos not in struct['plddt']:
        return res

    plddt = struct['plddt']
    seq = struct['sequence']
    coords_nz = struct['coords_nz']
    acidic_coords = struct['acidic_coords']
    acidic_positions = struct['acidic_positions']

    # pLDDT values
    res['plddt_site'] = plddt.get(pos)

    # Window average using configurable window length
    half_window = PLDDT_WINDOW_LENGTH // 2
    window_vals = [plddt[i] for i in range(pos - half_window, pos + half_window + 1) if i in plddt]
    res['plddt_window_avg'] = round(sum(window_vals) / len(window_vals), 2) if window_vals else None

    # Flanking amino acids and site itself
    res['aa_m2'] = seq.get(pos - 2)
    res['aa_m1'] = seq.get(pos - 1)
    res['aa_site'] = seq.get(pos)
    res['aa_p1'] = seq.get(pos + 1)
    res['aa_p2'] = seq.get(pos + 2)

    # Calculate distances using NZ of Lys to CG/CD of acidic residues
    if pos in coords_nz and acidic_positions:
        lys_nz_coord = coords_nz[pos]

        min_any, cls_any = float('inf'), None
        min_not_flank, cls_not_flank = float('inf'), None
        n_within = 0

        for ac_pos in acidic_positions:
            if ac_pos not in acidic_coords:
                continue
            d = calc_distance(lys_nz_coord, acidic_coords[ac_pos])
            rel_pos = ac_pos - pos

            # Count all acidic within threshold (including neighbors)
            if d <= DISTANCE_THRESHOLD:
                n_within += 1

            # Closest acidic anywhere (any position)
            if d < min_any:
                min_any, cls_any = d, ac_pos

            # Closest acidic NOT at -2 or +2 (but can be -1, +1, or elsewhere)
            if rel_pos != -2 and rel_pos != 2:
                if d < min_not_flank:
                    min_not_flank, cls_not_flank = d, ac_pos

        res['n_acidic_within_threshold'] = n_within

        if cls_any:
            res['dist_acidic_any'] = round(min_any, 2)
            res['closest_acidic_any'] = f"{seq.get(cls_any, '?')}{cls_any}"
        if cls_not_flank:
            res['dist_acidic_not_flank'] = round(min_not_flank, 2)
            res['closest_acidic_not_flank'] = f"{seq.get(cls_not_flank, '?')}{cls_not_flank}"

    # Flexibility classification based on window average pLDDT
    if res['plddt_window_avg'] is not None:
        if res['plddt_window_avg'] < PLDDT_THRESHOLD:
            res['flexible'] = 'Yes'
        else:
            res['structured'] = 'Yes'

    # Consensus motifs
    # Forward: ψKxE (hydrophobic at -1, acidic at +2)
    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    # Inverse: ExKψ (acidic at -2, hydrophobic at +1)
    res['inverse_consensus'] = 'Yes' if res['aa_m2'] in ACIDIC_RESIDUES and res['aa_p1'] in HYDROPHOBIC_RESIDUES else None

    # Acidic in -2/+2 specifically (not -1/+1)
    if res['aa_m2'] in ACIDIC_RESIDUES or res['aa_p2'] in ACIDIC_RESIDUES:
        res['acidic_in_flank'] = 'Yes'

    # Acidic in space within threshold (other than -2 or +2 specifically)
    if res['dist_acidic_not_flank'] is not None and res['dist_acidic_not_flank'] <= DISTANCE_THRESHOLD:
        res['acidic_in_space'] = 'Yes'

    # Any acidic within threshold (can be neighboring)
    if res['n_acidic_within_threshold'] > 0:
        res['any_acidic_within_threshold'] = 'Yes'

    # Category assignment (hierarchical: consensus > flankacidic > space > none)
    is_consensus = (res['forward_consensus'] == 'Yes' or res['inverse_consensus'] == 'Yes')
    is_flank_acidic = (res['acidic_in_flank'] == 'Yes')
    is_space_acidic = (res['acidic_in_space'] == 'Yes')
    has_any_acidic = (res['any_acidic_within_threshold'] == 'Yes')

    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"

    if is_consensus:
        suffix = "_consensus"
    elif is_flank_acidic:
        suffix = "_flankacidic"
    elif has_any_acidic:  # This covers space acidic and neighboring acidic
        suffix = "_space"
    else:
        suffix = "_none"

    res['category'] = prefix + suffix
    return res


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_file(input_file: str, output_file: str = None):
    """Process input Excel file and generate output with SUMO site analysis."""
    print(f"\n{'=' * 60}\nSUMO Site Analyzer - Core Output\n{'=' * 60}\n")
    print(f"Parameters:")
    print(f"  pLDDT threshold: {PLDDT_THRESHOLD}")
    print(f"  pLDDT window length: {PLDDT_WINDOW_LENGTH}")
    print(f"  Distance threshold: {DISTANCE_THRESHOLD} Å")
    print(f"  Score thresholds: {SCORE_VERY_HIGH_THRESHOLD}/{SCORE_HIGH_THRESHOLD}/{SCORE_MEDIUM_THRESHOLD}")
    print(f"  Sliding window length: {SLIDING_WINDOW_LENGTH}")
    print()

    df = pd.read_excel(input_file, sheet_name=0, header=1)
    df = df.iloc[:, :59]  # Keep columns A to BG

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

    # Build result DataFrame
    res_df = pd.DataFrame(results)

    # Calculate sliding window fractions
    # Create binary columns for calculations
    is_flexible = (res_df['flexible'] == 'Yes').astype(int)
    is_consensus = ((res_df['forward_consensus'] == 'Yes') | (res_df['inverse_consensus'] == 'Yes')).astype(int)
    has_any_acidic = (res_df['any_acidic_within_threshold'] == 'Yes').astype(int)

    # Calculate sliding window fractions (centered window)
    half_window = SLIDING_WINDOW_LENGTH // 2

    def calc_sliding_fraction(series):
        """Calculate sliding window fraction for a binary series."""
        n = len(series)
        fractions = []
        for i in range(n):
            start = max(0, i - half_window)
            end = min(n, i + half_window + 1)
            window = series.iloc[start:end]
            fractions.append(window.sum() / len(window) if len(window) > 0 else None)
        return fractions

    res_df['sliding_frac_flexible'] = calc_sliding_fraction(is_flexible)
    res_df['sliding_frac_consensus'] = calc_sliding_fraction(is_consensus)
    res_df['sliding_frac_any_acidic'] = calc_sliding_fraction(has_any_acidic)

    # Column mapping to final names
    cols_map = {
        'used_id': 'UniProt_ID_used',
        'mapped_pos': 'Mapped_position',
        'plddt_site': 'pLDDT_site',
        'plddt_window_avg': 'pLDDT_window_avg',
        'aa_m2': 'AA_-2',
        'aa_m1': 'AA_-1',
        'aa_site': 'AA_site',
        'aa_p1': 'AA_+1',
        'aa_p2': 'AA_+2',
        'dist_acidic_any': 'Dist_acidic_any_A',
        'closest_acidic_any': 'Closest_acidic_any',
        'dist_acidic_not_flank': 'Dist_acidic_notflank_A',
        'closest_acidic_not_flank': 'Closest_acidic_notflank',
        'n_acidic_within_threshold': 'N_acidic_within_threshold',
        'flexible': 'Flexible',
        'structured': 'Structured',
        'forward_consensus': 'Forward_consensus',
        'inverse_consensus': 'Inverse_consensus',
        'acidic_in_flank': 'Acidic_in_-2/+2',
        'acidic_in_space': 'Acidic_in_space',
        'any_acidic_within_threshold': 'Any_acidic_within_threshold',
        'category': 'Category',
        'sliding_frac_flexible': 'Sliding_frac_Flexible',
        'sliding_frac_consensus': 'Sliding_frac_Consensus',
        'sliding_frac_any_acidic': 'Sliding_frac_Any_acidic'
    }
    res_df = res_df.rename(columns=cols_map)

    # Order columns as requested
    ordered_cols = [
        'UniProt_ID_used',
        'Mapped_position',
        'pLDDT_site',
        'pLDDT_window_avg',
        'AA_-2',
        'AA_-1',
        'AA_site',
        'AA_+1',
        'AA_+2',
        'Dist_acidic_any_A',
        'Closest_acidic_any',
        'Dist_acidic_notflank_A',
        'Closest_acidic_notflank',
        'N_acidic_within_threshold',
        'Flexible',
        'Structured',
        'Forward_consensus',
        'Inverse_consensus',
        'Acidic_in_-2/+2',
        'Acidic_in_space',
        'Any_acidic_within_threshold',
        'Category',
        'Sliding_frac_Flexible',
        'Sliding_frac_Consensus',
        'Sliding_frac_Any_acidic'
    ]
    res_df = res_df[[c for c in ordered_cols if c in res_df.columns]]

    df = pd.concat([df, res_df], axis=1)

    # Summary
    total_sites = len(df)
    success_count = df['pLDDT_site'].notna().sum()
    failed_count = total_sites - success_count
    success_rate = (success_count / total_sites * 100) if total_sites > 0 else 0

    print(f"\n{'=' * 40}")
    print("DATA COMPLETENESS SUMMARY")
    print(f"{'=' * 40}")
    print(f"Total sites:           {total_sites}")
    print(f"Successfully resolved: {success_count} ({success_rate:.1f}%)")
    print(f"Failed to resolve:     {failed_count}")
    print(f"{'=' * 40}\n")

    # Output file
    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_core.xlsx'

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Sites', index=False)

    print(f"Output saved to: {output_file}")
    return df


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nConfigurable parameters (edit at top of script):")
        print(f"  PLDDT_THRESHOLD = {PLDDT_THRESHOLD}")
        print(f"  PLDDT_WINDOW_LENGTH = {PLDDT_WINDOW_LENGTH}")
        print(f"  DISTANCE_THRESHOLD = {DISTANCE_THRESHOLD}")
        print(f"  SLIDING_WINDOW_LENGTH = {SLIDING_WINDOW_LENGTH}")
        print(f"  SCORE_VERY_HIGH_THRESHOLD = {SCORE_VERY_HIGH_THRESHOLD}")
        print(f"  SCORE_HIGH_THRESHOLD = {SCORE_HIGH_THRESHOLD}")
        print(f"  SCORE_MEDIUM_THRESHOLD = {SCORE_MEDIUM_THRESHOLD}")
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
