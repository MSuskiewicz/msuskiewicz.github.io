#!/usr/bin/env python3
"""
SUMO Site Analyzer - Core Output (Exposed Acidic Focus)
========================================================
Parses Excel files with SUMO modification sites, fetches AlphaFold structures,
and generates core analysis data focusing on surface-exposed acidic residues.

Surface exposure is calculated using neighbor counting (approximates SASA).
Residues with fewer neighboring atoms are considered more solvent-exposed.

Usage:
    python sumo_site_core_exposed.py <input_excel_file> [output_file]
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

# Distance thresholds (separate for flexible and structured sites)
DISTANCE_THRESHOLD_FLEXIBLE = 8.0    # Angstroms for Cα-Cα distance in flexible sites
DISTANCE_THRESHOLD_STRUCTURED = 8.0  # Angstroms for NZ-CG/CD distance in structured sites

# Surface exposure parameters
EXPOSURE_NEIGHBOR_RADIUS = 10.0      # Radius (Å) to count neighboring Cβ atoms
EXPOSURE_THRESHOLD = 0.5             # Fraction: residue exposed if neighbors < (1-threshold) * max_neighbors
MAX_NEIGHBORS_BURIED = 24            # Approximate max neighbors for fully buried residue
# A residue is "exposed" if neighbor_count < MAX_NEIGHBORS_BURIED * (1 - EXPOSURE_THRESHOLD)
# With threshold=0.5 and max=24, exposed if < 12 neighbors

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


def calculate_exposure(coords_cb: dict, pos: int, all_cb_coords: list) -> bool:
    """
    Determine if a residue is surface-exposed based on neighbor counting.

    A residue is considered exposed if it has fewer neighboring Cβ atoms
    than expected for a buried residue.

    Args:
        coords_cb: Dict mapping residue number to Cβ coordinate (or Cα for Gly)
        pos: Residue position to check
        all_cb_coords: List of all Cβ coordinates in the structure

    Returns:
        True if residue is exposed, False if buried
    """
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

    # Exposed if fewer neighbors than threshold
    exposure_cutoff = MAX_NEIGHBORS_BURIED * (1 - EXPOSURE_THRESHOLD)
    return neighbor_count < exposure_cutoff


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
    """Parse mmCIF file to extract coordinates, pLDDT values, and exposure data.

    Extracts:
    - CA coordinates for all residues
    - CB coordinates for exposure calculation (CA for Gly)
    - NZ coordinates for lysines (for structured site distance calculations)
    - CG coordinates for Asp, CD coordinates for Glu (for structured site distances)
    - Acidic CA coordinates (for flexible site Cα-Cα distances)
    - pLDDT values from B-factor column
    """
    plddt, seq, coords_ca, coords_cb, coords_nz = {}, {}, {}, {}, {}
    acidic_coords_sidechain = {}  # CG for Asp, CD for Glu (for structured)
    acidic_coords_ca = {}  # CA for acidic residues (for flexible)
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
                    # For Gly, use CA as CB equivalent
                    if aa3 == 'GLY':
                        coords_cb[rn] = coord
                    # Store CA for acidic residues (for flexible site calculations)
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

    seq_str = ''.join(seq.get(i, 'X') for i in range(1, max(seq.keys(), default=0) + 1))
    acidic_positions = sorted(set(acidic_coords_sidechain.keys()) | set(acidic_coords_ca.keys()))

    # Calculate exposure for all residues
    all_cb_coords = list(coords_cb.values())
    exposure = {}
    for pos in coords_cb:
        exposure[pos] = calculate_exposure(coords_cb, pos, all_cb_coords)

    return {
        'plddt': plddt,
        'sequence': seq,
        'sequence_string': seq_str,
        'coords_ca': coords_ca,
        'coords_cb': coords_cb,
        'coords_nz': coords_nz,
        'acidic_coords_sidechain': acidic_coords_sidechain,  # CG/CD for structured
        'acidic_coords_ca': acidic_coords_ca,  # CA for flexible
        'acidic_positions': acidic_positions,
        'exposure': exposure  # Dict: position -> bool (True = exposed)
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
    """Analyze a SUMO site focusing on exposed acidic residues.

    Distance calculation depends on flexibility of BOTH residues:
    - If EITHER lysine OR acidic has pLDDT < threshold: use Cα-Cα distance
    - If BOTH lysine AND acidic have pLDDT >= threshold: use NZ-CG/CD distance

    Returns dict with all output columns including exposure data.
    """
    res = {
        'plddt_site': None,
        'plddt_window_avg': None,
        'aa_m2': None,
        'aa_m1': None,
        'aa_site': None,
        'aa_p1': None,
        'aa_p2': None,
        'exposed': None,  # Is the lysine site exposed?
        'dist_exposed_acidic': None,  # Distance to nearest exposed acidic
        'nearest_exposed_acidic': None,  # Position of nearest exposed acidic
        'close_exposed_acidic': None,  # Yes if exposed acidic within threshold
        'n_exposed_acidic_within': 0,  # Count of exposed acidics within threshold
        'flexible': None,
        'structured': None,
        'category': None  # Flexible_exposed, Flexible_noexposed, etc.
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

    # pLDDT values
    res['plddt_site'] = plddt.get(pos)
    lys_plddt = plddt.get(pos, 0)

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

    # Is lysine site exposed?
    res['exposed'] = 'Yes' if exposure.get(pos, False) else None

    # Determine flexibility based on window average
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
        min_exposed_dist = float('inf')
        nearest_exposed_pos = None
        n_exposed_within = 0

        for ac_pos in acidic_positions:
            # Only consider EXPOSED acidic residues
            if not exposure.get(ac_pos, False):
                continue

            ac_plddt = plddt.get(ac_pos, 0)

            # Determine which distance calculation to use for this pair
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

            # Count exposed acidics within threshold
            if d <= distance_threshold:
                n_exposed_within += 1

            # Track nearest exposed acidic
            if d < min_exposed_dist:
                min_exposed_dist = d
                nearest_exposed_pos = ac_pos

        res['n_exposed_acidic_within'] = n_exposed_within

        if nearest_exposed_pos is not None:
            res['dist_exposed_acidic'] = round(min_exposed_dist, 2)
            res['nearest_exposed_acidic'] = f"{seq.get(nearest_exposed_pos, '?')}{nearest_exposed_pos}"

        # Close exposed acidic if any within threshold
        if n_exposed_within > 0:
            res['close_exposed_acidic'] = 'Yes'

    # Category assignment based on exposed acidic presence
    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"
    suffix = "_exposed" if res['close_exposed_acidic'] == 'Yes' else "_noexposed"
    res['category'] = prefix + suffix

    return res


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_file(input_file: str, output_file: str = None):
    """Process input Excel file and generate output with SUMO site analysis (exposed focus)."""
    print(f"\n{'=' * 60}\nSUMO Site Analyzer - Core Output (Exposed Acidic Focus)\n{'=' * 60}\n")
    print(f"Parameters:")
    print(f"  pLDDT threshold: {PLDDT_THRESHOLD}")
    print(f"  pLDDT window length: {PLDDT_WINDOW_LENGTH}")
    print(f"  Distance threshold (flexible, Cα-Cα): {DISTANCE_THRESHOLD_FLEXIBLE} Å")
    print(f"  Distance threshold (structured, NZ-CG/CD): {DISTANCE_THRESHOLD_STRUCTURED} Å")
    print(f"  Exposure neighbor radius: {EXPOSURE_NEIGHBOR_RADIUS} Å")
    print(f"  Exposure threshold: {EXPOSURE_THRESHOLD} ({int(EXPOSURE_THRESHOLD*100)}%)")
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
    is_flexible = (res_df['flexible'] == 'Yes').astype(int)
    is_structured = (res_df['structured'] == 'Yes').astype(int)
    has_close_exposed = (res_df['close_exposed_acidic'] == 'Yes').astype(int)

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

    def calc_conditional_sliding_fraction(condition_series, value_series):
        """Calculate sliding window fraction of value_series among condition_series==True.

        For each position, looks at the surrounding window and calculates:
        (count of sites where BOTH condition AND value are true) / (count where condition is true)
        """
        n = len(condition_series)
        fractions = []
        for i in range(n):
            start = max(0, i - half_window)
            end = min(n, i + half_window + 1)
            cond_window = condition_series.iloc[start:end]
            val_window = value_series.iloc[start:end]
            n_cond = cond_window.sum()
            if n_cond > 0:
                n_both = (cond_window & val_window).sum()
                fractions.append(n_both / n_cond)
            else:
                fractions.append(None)
        return fractions

    res_df['sliding_frac_flexible'] = calc_sliding_fraction(is_flexible)
    res_df['sliding_frac_exposed_acidic'] = calc_sliding_fraction(has_close_exposed)

    # Conditional sliding fractions: fraction of Flex/Struct sites that have exposed acidic nearby
    res_df['sliding_frac_flex_with_acidic'] = calc_conditional_sliding_fraction(
        is_flexible == 1, has_close_exposed == 1)
    res_df['sliding_frac_struct_with_acidic'] = calc_conditional_sliding_fraction(
        is_structured == 1, has_close_exposed == 1)

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
        'exposed': 'Lys_exposed',
        'dist_exposed_acidic': 'Dist_exposed_acidic_A',
        'nearest_exposed_acidic': 'Nearest_exposed_acidic',
        'close_exposed_acidic': 'Close_exposed_acidic',
        'n_exposed_acidic_within': 'N_exposed_acidic_within',
        'flexible': 'Flexible',
        'structured': 'Structured',
        'category': 'Category',
        'sliding_frac_flexible': 'Sliding_frac_Flexible',
        'sliding_frac_exposed_acidic': 'Sliding_frac_Exposed_acidic',
        'sliding_frac_flex_with_acidic': 'Sliding_frac_Flex_with_acidic',
        'sliding_frac_struct_with_acidic': 'Sliding_frac_Struct_with_acidic'
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
        'Lys_exposed',
        'Dist_exposed_acidic_A',
        'Nearest_exposed_acidic',
        'Close_exposed_acidic',
        'N_exposed_acidic_within',
        'Flexible',
        'Structured',
        'Category',
        'Sliding_frac_Flexible',
        'Sliding_frac_Exposed_acidic',
        'Sliding_frac_Flex_with_acidic',
        'Sliding_frac_Struct_with_acidic'
    ]
    res_df = res_df[[c for c in ordered_cols if c in res_df.columns]]

    df = pd.concat([df, res_df], axis=1)

    # Summary
    total_sites = len(df)
    success_count = df['pLDDT_site'].notna().sum()
    failed_count = total_sites - success_count
    success_rate = (success_count / total_sites * 100) if total_sites > 0 else 0

    # Exposure summary
    exposed_count = (df['Lys_exposed'] == 'Yes').sum() if 'Lys_exposed' in df.columns else 0
    has_exposed_acidic = (df['Close_exposed_acidic'] == 'Yes').sum() if 'Close_exposed_acidic' in df.columns else 0

    print(f"\n{'=' * 40}")
    print("DATA COMPLETENESS SUMMARY")
    print(f"{'=' * 40}")
    print(f"Total sites:              {total_sites}")
    print(f"Successfully resolved:    {success_count} ({success_rate:.1f}%)")
    print(f"Failed to resolve:        {failed_count}")
    print(f"Exposed lysine sites:     {exposed_count}")
    print(f"Has close exposed acidic: {has_exposed_acidic}")
    print(f"{'=' * 40}\n")

    # Output file
    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_core_exposed.xlsx'

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
        print(f"  DISTANCE_THRESHOLD_FLEXIBLE = {DISTANCE_THRESHOLD_FLEXIBLE} (Cα-Cα)")
        print(f"  DISTANCE_THRESHOLD_STRUCTURED = {DISTANCE_THRESHOLD_STRUCTURED} (NZ-CG/CD)")
        print(f"  EXPOSURE_NEIGHBOR_RADIUS = {EXPOSURE_NEIGHBOR_RADIUS}")
        print(f"  EXPOSURE_THRESHOLD = {EXPOSURE_THRESHOLD}")
        print(f"  SLIDING_WINDOW_LENGTH = {SLIDING_WINDOW_LENGTH}")
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
