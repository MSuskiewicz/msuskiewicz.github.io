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
    """Parse mmCIF file to extract coordinates, pLDDT values, and secondary structure.

    Extracts:
    - CA coordinates for all residues (for Cα-Cα distance calculations)
    - CB coordinates for all residues (for exposure calculation via neighbor counting)
    - pLDDT values from B-factor column
    - Secondary structure from _struct_conf category (helix/strand/coil)
    - Acidic residue (D, E) coordinates
    - Hydrophobic residue coordinates
    """
    plddt, seq, coords_ca, coords_cb = {}, {}, {}, {}
    acidic_coords_ca = {}  # CA for acidic residues
    hydrophobic_coords_ca = {}  # CA for hydrophobic residues
    in_atom, headers, idx = False, [], {}

    # Secondary structure storage
    sec_struct = {}  # residue_num -> 'helix' | 'strand' | 'coil'

    # First pass: extract secondary structure from _struct_conf
    lines = content.split('\n')
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
            # Parse secondary structure entry
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
                    # Store CA for acidic residues
                    if aa3 in ('ASP', 'GLU'):
                        acidic_coords_ca[rn] = coord
                    # Store CA for hydrophobic residues
                    if aa in HYDROPHOBIC_RESIDUES:
                        hydrophobic_coords_ca[rn] = coord
                elif atom_name == 'CB':
                    coords_cb[rn] = coord
            except Exception:
                continue
        elif in_atom and line.startswith('#'):
            break

    seq_str = ''.join(seq.get(i, 'X') for i in range(1, max(seq.keys(), default=0) + 1))
    acidic_positions = sorted(acidic_coords_ca.keys())
    hydrophobic_positions = sorted(hydrophobic_coords_ca.keys())

    # For glycines (no CB), use CA as fallback for exposure calculation
    for rn in coords_ca:
        if rn not in coords_cb:
            coords_cb[rn] = coords_ca[rn]

    return {
        'plddt': plddt,
        'sequence': seq,
        'sequence_string': seq_str,
        'coords_ca': coords_ca,
        'coords_cb': coords_cb,  # For exposure calculation
        'acidic_coords_ca': acidic_coords_ca,
        'acidic_positions': acidic_positions,
        'hydrophobic_coords_ca': hydrophobic_coords_ca,
        'hydrophobic_positions': hydrophobic_positions,
        'secondary_structure': sec_struct  # helix/strand/coil
    }


# =============================================================================
# EXPOSURE CALCULATION
# =============================================================================

def calculate_exposure(struct: dict, pos: int) -> bool:
    """Determine if a residue is surface-exposed using neighbor counting.

    A residue is considered exposed if it has fewer than EXPOSURE_NEIGHBOR_THRESHOLD
    CB atoms within EXPOSURE_DISTANCE_THRESHOLD Angstroms.

    Returns True if exposed, False if buried.
    """
    coords_cb = struct.get('coords_cb', {})
    if pos not in coords_cb:
        return False  # Can't determine

    target_cb = coords_cb[pos]
    neighbor_count = 0

    for other_pos, other_cb in coords_cb.items():
        if other_pos == pos:
            continue
        d = calc_distance(target_cb, other_cb)
        if d <= EXPOSURE_DISTANCE_THRESHOLD:
            neighbor_count += 1

    return neighbor_count < EXPOSURE_NEIGHBOR_THRESHOLD


# =============================================================================
# RESOLUTION & ANALYSIS
# =============================================================================

def resolve_structure_and_position(uid: str, isoform: int | None, pos: int, cache: dict) -> tuple[dict | None, str, int, str | None, str | None]:
    """Resolve AlphaFold structure and map position if needed.

    Returns:
        (struct, used_id, final_pos, original_seq, af_seq)
        - original_seq and af_seq are returned when isoform mapping was used,
          to enable reverse mapping of acidic positions back to original coords
    """
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
        original_seq = None
        af_seq = None
        if isoform:
            iso_seq = get_uniprot_sequence(uid, isoform, cache)
            if iso_seq:
                mapped = map_position_by_residue(iso_seq, struct['sequence_string'], pos)
                if mapped:
                    final_pos = mapped
                    original_seq = iso_seq
                    af_seq = struct['sequence_string']
                else:
                    result = (None, uid, pos, None, None)
                    cache[cache_key] = result
                    return result
        result = (struct, uid, final_pos, original_seq, af_seq)
        cache[cache_key] = result
        return result

    original_seq = get_unisave_sequence(uid, cache, min_length=pos) or get_uniprot_sequence(uid, isoform, cache)
    if not original_seq or pos > len(original_seq):
        result = (None, uid, pos, None, None)
        cache[cache_key] = result
        return result

    target_acc = get_mapped_accession(uid, cache)
    if target_acc:
        struct = fetch_alphafold(target_acc, 1)
        if struct:
            mapped_pos = map_position_by_residue(original_seq, struct['sequence_string'], pos)
            if mapped_pos:
                result = (struct, target_acc, mapped_pos, original_seq, struct['sequence_string'])
                cache[cache_key] = result
                return result

    result = (None, uid, pos, None, None)
    cache[cache_key] = result
    return result


def analyze_site(struct: dict, pos: int, original_pos: int = None,
                 original_seq: str = None, af_seq: str = None) -> dict:
    """Analyze a SUMO site and compute all metrics.

    Uses Cα-Cα distance for all calculations.
    Distance thresholds: minimum 4.9Å, maximum 8.0Å.

    When original_seq and af_seq are provided, acidic positions are mapped back
    to original sequence coordinates for consistent reporting.

    Args:
        struct: Parsed AlphaFold structure
        pos: Position in AF model coordinates
        original_pos: Original position before mapping (for relative position calc)
        original_seq: Original isoform/input sequence (for reverse mapping)
        af_seq: AlphaFold model sequence (for reverse mapping)

    Returns dict with all output columns including:
    - Lists of acidic residues within threshold (with exposure markers)
    - Distances in Angstroms
    - Relative positions from lysine
    - Secondary structure
    """
    res = {
        'plddt_site': None,
        'plddt_window_avg': None,
        'aa_m2': None,
        'aa_m1': None,
        'aa_site': None,
        'aa_p1': None,
        'aa_p2': None,
        # Acidic analysis
        'acidics_within_threshold': None,      # List: "D123(exp),E456"
        'distances_angstroms': None,           # List: "5.2,6.8"
        'distances_positions': None,           # List: "-5(exp),+20"
        'any_acidic_within_threshold': None,
        'exposed_acidic_within_threshold': None,
        'acidic_in_pm2': None,                 # Acidic within +/-2 positions
        'exposed_acidic_in_pm2': None,         # Exposed acidic within +/-2 positions
        # Extended acidics (positions +3 to +10 and -3 to -10)
        'extended_acidics_downstream': 0,      # Count of acidics in +3 to +10
        'extended_acidics_upstream': 0,        # Count of acidics in -10 to -3
        'extended_acidics_total': 0,           # Sum of upstream + downstream
        # Hydrophobic analysis
        'hydrophobics_within_threshold': None, # List: "V123(exp),L456"
        'hydrophobic_distances_angstroms': None,  # List: "3.2,4.1"
        'hydrophobic_distances_positions': None,  # List: "-3(exp),+5"
        'any_hydrophobic_within_threshold': None,
        'exposed_hydrophobic_within_threshold': None,
        # Structure info
        'flexible': None,
        'structured': None,
        'secondary_structure': None,           # helix/strand/coil
        'forward_consensus': None,
        'inverse_consensus': None,
        'category': None
    }

    if not struct or pos not in struct['plddt']:
        return res

    # Use original position for relative calculations if provided
    site_pos_for_rel = original_pos if original_pos is not None else pos

    plddt = struct['plddt']
    seq = struct['sequence']
    coords_ca = struct['coords_ca']
    acidic_coords_ca = struct['acidic_coords_ca']
    acidic_positions = struct['acidic_positions']
    sec_struct = struct.get('secondary_structure', {})

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

    # Determine flexibility based on window average
    if res['plddt_window_avg'] is not None:
        if res['plddt_window_avg'] < PLDDT_THRESHOLD:
            res['flexible'] = 'Yes'
        else:
            res['structured'] = 'Yes'

    # Secondary structure from AlphaFold model
    ss = sec_struct.get(pos)
    if ss:
        res['secondary_structure'] = ss
    else:
        res['secondary_structure'] = 'coil'  # Default if not defined

    # Calculate Cα-Cα distances to all acidic residues
    lys_ca = coords_ca.get(pos)

    if acidic_positions and lys_ca:
        # Collect all acidics within the threshold range (min to max)
        # Store: (acidic_pos_af, acidic_pos_orig, distance, is_exposed)
        acidics_info = []

        for ac_pos in acidic_positions:
            if ac_pos not in acidic_coords_ca:
                continue

            d = calc_distance(lys_ca, acidic_coords_ca[ac_pos])

            # Check if within threshold range (min <= d <= max)
            if DISTANCE_THRESHOLD_MIN <= d <= DISTANCE_THRESHOLD_MAX:
                is_exposed = calculate_exposure(struct, ac_pos)

                # Map acidic position back to original coordinates if needed
                if original_seq and af_seq:
                    # Reverse map: AF model position -> original isoform position
                    mapped_ac_pos = map_position_by_residue(af_seq, original_seq, ac_pos)
                    if mapped_ac_pos is None:
                        # This acidic doesn't exist in the original isoform - skip it
                        continue
                    ac_pos_for_report = mapped_ac_pos
                else:
                    ac_pos_for_report = ac_pos

                acidics_info.append((ac_pos, ac_pos_for_report, d, is_exposed))

        # Sort by distance
        acidics_info.sort(key=lambda x: x[2])

        if acidics_info:
            # Build lists for output
            acidic_labels = []
            distance_labels = []
            position_labels = []

            has_exposed = False
            has_in_pm2 = False           # Any acidic within +/-2 positions
            has_exposed_in_pm2 = False   # Exposed acidic within +/-2 positions

            for ac_pos_af, ac_pos_report, dist, is_exposed in acidics_info:
                aa = seq.get(ac_pos_af, '?')
                # Use reported position for relative calculation
                rel_pos = ac_pos_report - site_pos_for_rel
                rel_pos_str = f"+{rel_pos}" if rel_pos > 0 else str(rel_pos)

                # Check if within +/-2 positions
                if -2 <= rel_pos <= 2:
                    has_in_pm2 = True
                    if is_exposed:
                        has_exposed_in_pm2 = True

                if is_exposed:
                    acidic_labels.append(f"{aa}{ac_pos_report}(exp)")
                    distance_labels.append(f"{dist:.1f}(exp)")
                    position_labels.append(f"{rel_pos_str}(exp)")
                    has_exposed = True
                else:
                    acidic_labels.append(f"{aa}{ac_pos_report}")
                    distance_labels.append(f"{dist:.1f}")
                    position_labels.append(rel_pos_str)

            res['acidics_within_threshold'] = ','.join(acidic_labels)
            res['distances_angstroms'] = ','.join(distance_labels)
            res['distances_positions'] = ','.join(position_labels)
            res['any_acidic_within_threshold'] = 'Yes'

            if has_exposed:
                res['exposed_acidic_within_threshold'] = 'Yes'
            if has_in_pm2:
                res['acidic_in_pm2'] = 'Yes'
            if has_exposed_in_pm2:
                res['exposed_acidic_in_pm2'] = 'Yes'

    # Calculate extended acidics (based on sequence position, not distance threshold)
    # Downstream: positions +3 to +10 relative to lysine
    # Upstream: positions -10 to -3 relative to lysine
    extended_downstream = 0
    extended_upstream = 0

    for offset in range(3, 11):  # +3 to +10
        check_pos = pos + offset
        aa_at_pos = seq.get(check_pos)
        if aa_at_pos in ACIDIC_RESIDUES:
            extended_downstream += 1

    for offset in range(-10, -2):  # -10 to -3
        check_pos = pos + offset
        aa_at_pos = seq.get(check_pos)
        if aa_at_pos in ACIDIC_RESIDUES:
            extended_upstream += 1

    res['extended_acidics_downstream'] = extended_downstream
    res['extended_acidics_upstream'] = extended_upstream
    res['extended_acidics_total'] = extended_downstream + extended_upstream

    # Calculate Cα-Cα distances to all hydrophobic residues
    hydrophobic_coords_ca = struct.get('hydrophobic_coords_ca', {})
    hydrophobic_positions = struct.get('hydrophobic_positions', [])

    if hydrophobic_positions and lys_ca:
        # Collect all hydrophobics within the threshold range
        hydrophobics_info = []

        for hp_pos in hydrophobic_positions:
            if hp_pos not in hydrophobic_coords_ca:
                continue

            d = calc_distance(lys_ca, hydrophobic_coords_ca[hp_pos])

            # Check if within hydrophobic threshold range
            if HYDROPHOBIC_DISTANCE_MIN <= d <= HYDROPHOBIC_DISTANCE_MAX:
                is_exposed = calculate_exposure(struct, hp_pos)

                # Map position back to original coordinates if needed
                if original_seq and af_seq:
                    mapped_hp_pos = map_position_by_residue(af_seq, original_seq, hp_pos)
                    if mapped_hp_pos is None:
                        continue
                    hp_pos_for_report = mapped_hp_pos
                else:
                    hp_pos_for_report = hp_pos

                hydrophobics_info.append((hp_pos, hp_pos_for_report, d, is_exposed))

        # Sort by distance
        hydrophobics_info.sort(key=lambda x: x[2])

        if hydrophobics_info:
            hp_labels = []
            hp_distance_labels = []
            hp_position_labels = []
            has_exposed_hp = False

            for hp_pos_af, hp_pos_report, dist, is_exposed in hydrophobics_info:
                aa = seq.get(hp_pos_af, '?')
                # Calculate relative position from lysine
                rel_pos = hp_pos_report - site_pos_for_rel
                rel_pos_str = f"+{rel_pos}" if rel_pos > 0 else str(rel_pos)

                if is_exposed:
                    hp_labels.append(f"{aa}{hp_pos_report}(exp)")
                    hp_distance_labels.append(f"{dist:.1f}(exp)")
                    hp_position_labels.append(f"{rel_pos_str}(exp)")
                    has_exposed_hp = True
                else:
                    hp_labels.append(f"{aa}{hp_pos_report}")
                    hp_distance_labels.append(f"{dist:.1f}")
                    hp_position_labels.append(rel_pos_str)

            res['hydrophobics_within_threshold'] = ','.join(hp_labels)
            res['hydrophobic_distances_angstroms'] = ','.join(hp_distance_labels)
            res['hydrophobic_distances_positions'] = ','.join(hp_position_labels)
            res['any_hydrophobic_within_threshold'] = 'Yes'

            if has_exposed_hp:
                res['exposed_hydrophobic_within_threshold'] = 'Yes'

    # Consensus motifs
    # Forward: ψKxE (hydrophobic at -1, acidic at +2)
    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    # Inverse: ExKψ (acidic at -2, hydrophobic at +1)
    res['inverse_consensus'] = 'Yes' if res['aa_m2'] in ACIDIC_RESIDUES and res['aa_p1'] in HYDROPHOBIC_RESIDUES else None

    # Category assignment
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
    print(f"  Acidic distance threshold (Cα-Cα): {DISTANCE_THRESHOLD_MIN} - {DISTANCE_THRESHOLD_MAX} Å")
    print(f"  Hydrophobic distance threshold (Cα-Cα): {HYDROPHOBIC_DISTANCE_MIN} - {HYDROPHOBIC_DISTANCE_MAX} Å")
    print(f"  Exposure parameters: {EXPOSURE_DISTANCE_THRESHOLD} Å radius, <{EXPOSURE_NEIGHBOR_THRESHOLD} neighbors = exposed")
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

        struct, used_id, final_pos, original_seq, af_seq = resolve_structure_and_position(uid, iso, pos, cache)
        # Pass original position and sequences for consistent coordinate reporting
        analysis = analyze_site(struct, final_pos, original_pos=pos,
                                original_seq=original_seq, af_seq=af_seq)
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
    is_structured = (res_df['structured'] == 'Yes').astype(int)
    is_consensus = ((res_df['forward_consensus'] == 'Yes') | (res_df['inverse_consensus'] == 'Yes')).astype(int)
    has_any_acidic = (res_df['any_acidic_within_threshold'] == 'Yes').astype(int)
    has_exposed_acidic = (res_df['exposed_acidic_within_threshold'] == 'Yes').astype(int)

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

    def calc_conditional_sliding_fraction(condition_series, value_series):
        """Calculate sliding window fraction of value_series among condition_series==True."""
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
    res_df['sliding_frac_consensus'] = calc_sliding_fraction(is_consensus)
    res_df['sliding_frac_any_acidic'] = calc_sliding_fraction(has_any_acidic)
    res_df['sliding_frac_exposed_acidic'] = calc_sliding_fraction(has_exposed_acidic)

    # Fraction of Flexible sites that have acidic (among Flexible sites only)
    res_df['sliding_frac_flex_with_acidic'] = calc_conditional_sliding_fraction(
        is_flexible == 1, has_any_acidic == 1)
    # Fraction of Structured sites that have acidic (among Structured sites only)
    res_df['sliding_frac_struct_with_acidic'] = calc_conditional_sliding_fraction(
        is_structured == 1, has_any_acidic == 1)

    # Sliding window averages for extended acidics
    def calc_sliding_mean(series):
        """Calculate sliding window mean for a numeric series."""
        n = len(series)
        means = []
        for i in range(n):
            start = max(0, i - half_window)
            end = min(n, i + half_window + 1)
            window = series.iloc[start:end]
            means.append(window.mean() if len(window) > 0 else None)
        return means

    res_df['sliding_avg_ext_acidics_downstream'] = calc_sliding_mean(res_df['extended_acidics_downstream'])
    res_df['sliding_avg_ext_acidics_upstream'] = calc_sliding_mean(res_df['extended_acidics_upstream'])
    res_df['sliding_avg_ext_acidics_total'] = calc_sliding_mean(res_df['extended_acidics_total'])

    # Sliding window fraction for hydrophobic
    has_any_hydrophobic = (res_df['any_hydrophobic_within_threshold'] == 'Yes').astype(int)
    res_df['sliding_frac_any_hydrophobic'] = calc_sliding_fraction(has_any_hydrophobic)

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
        # Acidic columns
        'acidics_within_threshold': 'Acidics_within_threshold',
        'distances_angstroms': 'Distances_Angstroms',
        'distances_positions': 'Distances_positions',
        'any_acidic_within_threshold': 'Any_acidic_within_threshold',
        'exposed_acidic_within_threshold': 'Exposed_acidic_within_threshold',
        'acidic_in_pm2': 'Acidic_in_pm2',
        'exposed_acidic_in_pm2': 'Exposed_acidic_in_pm2',
        # Extended acidics columns
        'extended_acidics_downstream': 'Extended_acidics_downstream',
        'extended_acidics_upstream': 'Extended_acidics_upstream',
        'extended_acidics_total': 'Extended_acidics_total',
        # Hydrophobic columns
        'hydrophobics_within_threshold': 'Hydrophobics_within_threshold',
        'hydrophobic_distances_angstroms': 'Hydrophobic_distances_Angstroms',
        'hydrophobic_distances_positions': 'Hydrophobic_distances_positions',
        'any_hydrophobic_within_threshold': 'Any_hydrophobic_within_threshold',
        'exposed_hydrophobic_within_threshold': 'Exposed_hydrophobic_within_threshold',
        # Structure columns
        'flexible': 'Flexible',
        'structured': 'Structured',
        'secondary_structure': 'Secondary_structure',
        'forward_consensus': 'Forward_consensus',
        'inverse_consensus': 'Inverse_consensus',
        'category': 'Category',
        # Sliding window columns
        'sliding_frac_flexible': 'Sliding_frac_Flexible',
        'sliding_frac_consensus': 'Sliding_frac_Consensus',
        'sliding_frac_any_acidic': 'Sliding_frac_Any_acidic',
        'sliding_frac_exposed_acidic': 'Sliding_frac_Exposed_acidic',
        'sliding_frac_flex_with_acidic': 'Sliding_frac_Flex_with_acidic',
        'sliding_frac_struct_with_acidic': 'Sliding_frac_Struct_with_acidic',
        'sliding_avg_ext_acidics_downstream': 'Sliding_avg_Ext_acidics_downstream',
        'sliding_avg_ext_acidics_upstream': 'Sliding_avg_Ext_acidics_upstream',
        'sliding_avg_ext_acidics_total': 'Sliding_avg_Ext_acidics_total',
        'sliding_frac_any_hydrophobic': 'Sliding_frac_Any_hydrophobic'
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
        # Acidic columns
        'Acidics_within_threshold',
        'Distances_Angstroms',
        'Distances_positions',
        'Any_acidic_within_threshold',
        'Exposed_acidic_within_threshold',
        'Acidic_in_pm2',
        'Exposed_acidic_in_pm2',
        # Extended acidics columns
        'Extended_acidics_downstream',
        'Extended_acidics_upstream',
        'Extended_acidics_total',
        # Hydrophobic columns
        'Hydrophobics_within_threshold',
        'Hydrophobic_distances_Angstroms',
        'Hydrophobic_distances_positions',
        'Any_hydrophobic_within_threshold',
        'Exposed_hydrophobic_within_threshold',
        # Structure columns
        'Flexible',
        'Structured',
        'Secondary_structure',
        'Forward_consensus',
        'Inverse_consensus',
        'Category',
        # Sliding window columns
        'Sliding_frac_Flexible',
        'Sliding_frac_Consensus',
        'Sliding_frac_Any_acidic',
        'Sliding_frac_Exposed_acidic',
        'Sliding_frac_Flex_with_acidic',
        'Sliding_frac_Struct_with_acidic',
        'Sliding_avg_Ext_acidics_downstream',
        'Sliding_avg_Ext_acidics_upstream',
        'Sliding_avg_Ext_acidics_total',
        'Sliding_frac_Any_hydrophobic'
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
        print(f"  DISTANCE_THRESHOLD_MIN = {DISTANCE_THRESHOLD_MIN} Å (Cα-Cα, acidic)")
        print(f"  DISTANCE_THRESHOLD_MAX = {DISTANCE_THRESHOLD_MAX} Å (Cα-Cα, acidic)")
        print(f"  HYDROPHOBIC_DISTANCE_MIN = {HYDROPHOBIC_DISTANCE_MIN} Å (Cα-Cα)")
        print(f"  HYDROPHOBIC_DISTANCE_MAX = {HYDROPHOBIC_DISTANCE_MAX} Å (Cα-Cα)")
        print(f"  EXPOSURE_DISTANCE_THRESHOLD = {EXPOSURE_DISTANCE_THRESHOLD} Å")
        print(f"  EXPOSURE_NEIGHBOR_THRESHOLD = {EXPOSURE_NEIGHBOR_THRESHOLD}")
        print(f"  SLIDING_WINDOW_LENGTH = {SLIDING_WINDOW_LENGTH}")
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)


if __name__ == "__main__":
    main()
