#!/usr/bin/env python3
"""
SUMO Site Analyzer - Core Output (Simplest Version)
====================================================
Minimal version without threading, session pooling, or sliding windows.
Processes Excel files with SUMO modification sites, fetches AlphaFold structures,
and generates core analysis data.

Usage:
    python sumo_site_core_simplest.py <input_excel_file> [output_file]
"""

import sys
import os
import re
import math
import requests
import pandas as pd

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

PLDDT_THRESHOLD = 65.0
PLDDT_WINDOW_LENGTH = 11
DISTANCE_THRESHOLD_MIN = 4.9
DISTANCE_THRESHOLD_MAX = 8.0
HYDROPHOBIC_DISTANCE_MIN = 3.0
HYDROPHOBIC_DISTANCE_MAX = 4.5
EXPOSURE_DISTANCE_THRESHOLD = 10.0
EXPOSURE_NEIGHBOR_THRESHOLD = 18

HYDROPHOBIC_RESIDUES = {'A', 'V', 'L', 'I', 'M', 'F', 'C', 'P', 'Y'}
ACIDIC_RESIDUES = {'D', 'E'}

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def safe_get(url, timeout=30):
    """Simple HTTP GET with error handling."""
    try:
        r = requests.get(url, timeout=timeout)
        return r if r.status_code == 200 else None
    except:
        return None


def calc_distance(c1, c2):
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def parse_uniprot_id(protein_id):
    """Extract UniProt accession and isoform from protein ID."""
    match = re.search(r'([A-Z][0-9][A-Z0-9]{3,9}[0-9])(?:-(\d+))?', protein_id)
    if match:
        return match.group(1), int(match.group(2)) if match.group(2) else None
    return protein_id, None


def fetch_alphafold(uid):
    """Fetch AlphaFold structure."""
    r = safe_get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}")
    if not r:
        return None
    try:
        data = r.json()
        if not data:
            return None
        cif_url = data[0].get('cifUrl')
        if not cif_url:
            return None
        cif_r = safe_get(cif_url)
        return parse_cif(cif_r.text) if cif_r else None
    except:
        return None


def parse_cif(content):
    """Parse mmCIF to extract coordinates and pLDDT."""
    plddt, seq, coords_ca, coords_cb = {}, {}, {}, {}
    acidic_coords_ca, hydrophobic_coords_ca = {}, {}
    in_atom, headers, idx = False, [], {}

    lines = content.split('\n')

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
            except:
                continue
        elif in_atom and line.startswith('#'):
            break

    # For glycines, use CA as CB
    for rn in coords_ca:
        if rn not in coords_cb:
            coords_cb[rn] = coords_ca[rn]

    return {
        'plddt': plddt, 'sequence': seq, 'coords_ca': coords_ca, 'coords_cb': coords_cb,
        'acidic_coords_ca': acidic_coords_ca, 'acidic_positions': sorted(acidic_coords_ca.keys()),
        'hydrophobic_coords_ca': hydrophobic_coords_ca, 'hydrophobic_positions': sorted(hydrophobic_coords_ca.keys())
    }


def calculate_exposure(struct, pos):
    """Check if residue is surface-exposed."""
    coords_cb = struct.get('coords_cb', {})
    if pos not in coords_cb:
        return False
    target_cb = coords_cb[pos]
    neighbor_count = sum(1 for other_pos, other_cb in coords_cb.items()
                         if other_pos != pos and calc_distance(target_cb, other_cb) <= EXPOSURE_DISTANCE_THRESHOLD)
    return neighbor_count < EXPOSURE_NEIGHBOR_THRESHOLD


def analyze_site(struct, pos):
    """Analyze a SUMO site."""
    res = {
        'plddt_site': None, 'plddt_window_avg': None,
        'aa_m8': None, 'aa_m7': None, 'aa_m6': None, 'aa_m5': None, 'aa_m4': None, 'aa_m3': None,
        'aa_m2': None, 'aa_m1': None, 'aa_site': None,
        'aa_p1': None, 'aa_p2': None, 'aa_p3': None, 'aa_p4': None, 'aa_p5': None, 'aa_p6': None, 'aa_p7': None, 'aa_p8': None,
        'acidics_within_threshold': None, 'distances_angstroms': None, 'distances_positions': None,
        'any_acidic_within_threshold': None, 'exposed_acidic_within_threshold': None,
        'acidic_in_pm2': None, 'exposed_acidic_in_pm2': None,
        'hydrophobics_within_threshold': None, 'hydrophobic_distances_angstroms': None, 'hydrophobic_distances_positions': None,
        'any_hydrophobic_within_threshold': None, 'exposed_hydrophobic_within_threshold': None,
        'flexible': None, 'structured': None,
        'forward_consensus': None, 'inverse_consensus': None, 'category': None
    }

    if not struct or pos not in struct['plddt']:
        return res

    plddt, seq, coords_ca = struct['plddt'], struct['sequence'], struct['coords_ca']
    acidic_coords_ca, acidic_positions = struct['acidic_coords_ca'], struct['acidic_positions']

    # pLDDT
    res['plddt_site'] = plddt.get(pos)
    half_window = PLDDT_WINDOW_LENGTH // 2
    window_vals = [plddt[i] for i in range(pos - half_window, pos + half_window + 1) if i in plddt]
    res['plddt_window_avg'] = round(sum(window_vals) / len(window_vals), 2) if window_vals else None

    # Amino acids from -8 to +8
    for offset in range(-8, 9):
        key = f'aa_m{abs(offset)}' if offset < 0 else ('aa_site' if offset == 0 else f'aa_p{offset}')
        res[key] = seq.get(pos + offset)

    # Flexibility
    if res['plddt_window_avg'] is not None:
        res['flexible'] = 'Yes' if res['plddt_window_avg'] < PLDDT_THRESHOLD else None
        res['structured'] = 'Yes' if res['plddt_window_avg'] >= PLDDT_THRESHOLD else None

    # Acidic residues within threshold
    lys_ca = coords_ca.get(pos)
    if acidic_positions and lys_ca:
        acidics_info = []
        for ac_pos in acidic_positions:
            if ac_pos not in acidic_coords_ca:
                continue
            d = calc_distance(lys_ca, acidic_coords_ca[ac_pos])
            if DISTANCE_THRESHOLD_MIN <= d <= DISTANCE_THRESHOLD_MAX:
                is_exposed = calculate_exposure(struct, ac_pos)
                acidics_info.append((ac_pos, d, is_exposed))

        acidics_info.sort(key=lambda x: x[1])

        if acidics_info:
            acidic_labels, distance_labels, position_labels = [], [], []
            has_exposed, has_in_pm2, has_exposed_in_pm2 = False, False, False

            for ac_pos, dist, is_exposed in acidics_info:
                aa = seq.get(ac_pos, '?')
                rel_pos = ac_pos - pos
                rel_str = f"+{rel_pos}" if rel_pos > 0 else str(rel_pos)

                if -2 <= rel_pos <= 2:
                    has_in_pm2 = True
                    if is_exposed:
                        has_exposed_in_pm2 = True

                if is_exposed:
                    acidic_labels.append(f"{aa}{ac_pos}(exp)")
                    distance_labels.append(f"{dist:.1f}(exp)")
                    position_labels.append(f"{rel_str}(exp)")
                    has_exposed = True
                else:
                    acidic_labels.append(f"{aa}{ac_pos}")
                    distance_labels.append(f"{dist:.1f}")
                    position_labels.append(rel_str)

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

    # Hydrophobic residues
    hydrophobic_coords_ca = struct.get('hydrophobic_coords_ca', {})
    hydrophobic_positions = struct.get('hydrophobic_positions', [])
    if hydrophobic_positions and lys_ca:
        hp_info = []
        for hp_pos in hydrophobic_positions:
            if hp_pos not in hydrophobic_coords_ca:
                continue
            d = calc_distance(lys_ca, hydrophobic_coords_ca[hp_pos])
            if HYDROPHOBIC_DISTANCE_MIN <= d <= HYDROPHOBIC_DISTANCE_MAX:
                is_exposed = calculate_exposure(struct, hp_pos)
                hp_info.append((hp_pos, d, is_exposed))

        hp_info.sort(key=lambda x: x[1])

        if hp_info:
            hp_labels, hp_dist_labels, hp_pos_labels = [], [], []
            has_exposed_hp = False

            for hp_pos, dist, is_exposed in hp_info:
                aa = seq.get(hp_pos, '?')
                rel_pos = hp_pos - pos
                rel_str = f"+{rel_pos}" if rel_pos > 0 else str(rel_pos)

                if is_exposed:
                    hp_labels.append(f"{aa}{hp_pos}(exp)")
                    hp_dist_labels.append(f"{dist:.1f}(exp)")
                    hp_pos_labels.append(f"{rel_str}(exp)")
                    has_exposed_hp = True
                else:
                    hp_labels.append(f"{aa}{hp_pos}")
                    hp_dist_labels.append(f"{dist:.1f}")
                    hp_pos_labels.append(rel_str)

            res['hydrophobics_within_threshold'] = ','.join(hp_labels)
            res['hydrophobic_distances_angstroms'] = ','.join(hp_dist_labels)
            res['hydrophobic_distances_positions'] = ','.join(hp_pos_labels)
            res['any_hydrophobic_within_threshold'] = 'Yes'
            if has_exposed_hp:
                res['exposed_hydrophobic_within_threshold'] = 'Yes'

    # Consensus motifs
    res['forward_consensus'] = 'Yes' if res['aa_m1'] in HYDROPHOBIC_RESIDUES and res['aa_p2'] in ACIDIC_RESIDUES else None
    res['inverse_consensus'] = 'Yes' if res['aa_m2'] in ACIDIC_RESIDUES and res['aa_p1'] in HYDROPHOBIC_RESIDUES else None

    # Category - 5-level hierarchy
    is_consensus = res['forward_consensus'] == 'Yes' or res['inverse_consensus'] == 'Yes'
    has_exposed_pm2 = res['exposed_acidic_in_pm2'] == 'Yes'
    has_exposed = res['exposed_acidic_within_threshold'] == 'Yes'
    has_any = res['any_acidic_within_threshold'] == 'Yes'

    prefix = "Flexible" if res['flexible'] == 'Yes' else "Structured"
    if is_consensus:
        suffix = "_consensus"
    elif has_exposed_pm2:
        suffix = "_exposed_acidic_pm2"
    elif has_exposed:
        suffix = "_exposed_acidic"
    elif has_any:
        suffix = "_buried_acidic"
    else:
        suffix = "_no_acidic"
    res['category'] = prefix + suffix

    return res


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_file(input_file, output_file=None):
    """Process input Excel file."""
    print(f"Processing {input_file}...")

    df = pd.read_excel(input_file, sheet_name=0, header=1)
    df = df.iloc[:, :59]

    protein_col = next((c for c in df.columns if str(c).lower() == 'protein'), 'Protein')
    position_col = next((c for c in df.columns if str(c).lower() == 'position'), 'Position')

    results = []
    total = len(df)

    for idx, row in df.iterrows():
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{total}")

        uid, iso = parse_uniprot_id(str(row[protein_col]).strip())
        try:
            pos = int(row[position_col])
        except:
            results.append({'used_id': uid, 'mapped_pos': None})
            continue

        struct = fetch_alphafold(uid)
        analysis = analyze_site(struct, pos) if struct else {}
        analysis['used_id'] = uid
        analysis['mapped_pos'] = None
        results.append(analysis)

    res_df = pd.DataFrame(results)

    # Rename columns
    cols_map = {
        'used_id': 'UniProt_ID_used', 'mapped_pos': 'Mapped_position',
        'plddt_site': 'pLDDT_site', 'plddt_window_avg': 'pLDDT_window_avg',
        'aa_m8': 'AA_-8', 'aa_m7': 'AA_-7', 'aa_m6': 'AA_-6', 'aa_m5': 'AA_-5',
        'aa_m4': 'AA_-4', 'aa_m3': 'AA_-3', 'aa_m2': 'AA_-2', 'aa_m1': 'AA_-1',
        'aa_site': 'AA_site', 'aa_p1': 'AA_+1', 'aa_p2': 'AA_+2', 'aa_p3': 'AA_+3',
        'aa_p4': 'AA_+4', 'aa_p5': 'AA_+5', 'aa_p6': 'AA_+6', 'aa_p7': 'AA_+7', 'aa_p8': 'AA_+8',
        'acidics_within_threshold': 'Acidics_within_threshold',
        'distances_angstroms': 'Distances_Angstroms', 'distances_positions': 'Distances_positions',
        'any_acidic_within_threshold': 'Any_acidic_within_threshold',
        'exposed_acidic_within_threshold': 'Exposed_acidic_within_threshold',
        'acidic_in_pm2': 'Acidic_in_pm2', 'exposed_acidic_in_pm2': 'Exposed_acidic_in_pm2',
        'hydrophobics_within_threshold': 'Hydrophobics_within_threshold',
        'hydrophobic_distances_angstroms': 'Hydrophobic_distances_Angstroms',
        'hydrophobic_distances_positions': 'Hydrophobic_distances_positions',
        'any_hydrophobic_within_threshold': 'Any_hydrophobic_within_threshold',
        'exposed_hydrophobic_within_threshold': 'Exposed_hydrophobic_within_threshold',
        'flexible': 'Flexible', 'structured': 'Structured',
        'forward_consensus': 'Forward_consensus', 'inverse_consensus': 'Inverse_consensus',
        'category': 'Category'
    }
    res_df = res_df.rename(columns=cols_map)

    df = pd.concat([df, res_df], axis=1)

    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_core_simplest.xlsx'

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Sites', index=False)

    print(f"Output saved to: {output_file}")
    return df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
