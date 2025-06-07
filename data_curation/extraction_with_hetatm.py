import os
import logging
import requests
import time
import csv
import re
import statistics
from typing import Tuple, Optional, Dict, Any
from bs4 import BeautifulSoup
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import pandas as pd
import json

# ------------------------------------------------
# Constants / Configuration
# ------------------------------------------------
BASE_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json="
BASE_DOWNLOAD_URL = "https://files.rcsb.org/download/"
BASE_METADATA_URL = "https://data.rcsb.org/rest/v1/core/entry/"
BASE_STRUCTURE_URL = "https://www.rcsb.org/structure/"

NOT_ALLOWED_PDB_IDS = ["1KB2", "1KB4", "1KB6", "1YNW"]

OUTPUT_DIR = "pdb_files"
CSV_FILENAME = "vdr_structures_enhanced.csv"

SPECIES_LIST = [
    "Homo sapiens",
    "Rattus norvegicus", 
    "Danio rerio",
    "Petromyzon marinus",
]

RESOLUTION_BINS = {
    "1-2 Å": lambda r: 1.0 <= r < 2.0,
    "2-2.5 Å": lambda r: 2.0 <= r < 2.5,
    "2.5+ Å": lambda r: r >= 2.5
}

REQUEST_DELAY = 0.5
MAX_RETRIES = 3
RETRY_DELAY = 5

VDR_IDENTIFIERS = [
    "vitamin d receptor",
    "vdr",
    "vitamin d3 receptor",
    "calcitriol receptor",
    "vitamin nuclear receptor",
    "vitamin d nuclear receptor"
]

CRYSTAL_SYSTEMS = {
    'P 1': 'Triclinic',
    'P 2': 'Monoclinic',
    'P 21': 'Monoclinic',
    'C 2': 'Monoclinic',
    'P 1 21 1': 'Monoclinic',
    'C 1 2 1': 'Monoclinic',
    'P 2 2 2': 'Orthorhombic',
    'P 2 2 21': 'Orthorhombic',
    'P 21 21 2': 'Orthorhombic',
    'P 21 21 21': 'Orthorhombic',
    'C 2 2 21': 'Orthorhombic',
    'C 2 2 2': 'Orthorhombic',
    'F 2 2 2': 'Orthorhombic',
    'I 2 2 2': 'Orthorhombic',
    'I 21 21 21': 'Orthorhombic',
    'P 4': 'Tetragonal',
    'P 41': 'Tetragonal',
    'P 42': 'Tetragonal',
    'P 43': 'Tetragonal',
    'I 4': 'Tetragonal',
    'I 41': 'Tetragonal',
    'P 4 2 2': 'Tetragonal',
    'P 4 21 2': 'Tetragonal',
    'P 41 2 2': 'Tetragonal',
    'P 41 21 2': 'Tetragonal',
    'P 42 2 2': 'Tetragonal',
    'P 42 21 2': 'Tetragonal',
    'P 43 2 2': 'Tetragonal',
    'P 43 21 2': 'Tetragonal',
    'I 4 2 2': 'Tetragonal',
    'I 41 2 2': 'Tetragonal',
    'P 3': 'Trigonal',
    'P 31': 'Trigonal',
    'P 32': 'Trigonal',
    'R 3': 'Trigonal',
    'P 3 1 2': 'Trigonal',
    'P 3 2 1': 'Trigonal',
    'P 31 1 2': 'Trigonal',
    'P 31 2 1': 'Trigonal',
    'P 32 1 2': 'Trigonal',
    'P 32 2 1': 'Trigonal',
    'R 3 2': 'Trigonal',
    'P 6': 'Hexagonal',
    'P 61': 'Hexagonal',
    'P 65': 'Hexagonal',
    'P 62': 'Hexagonal',
    'P 64': 'Hexagonal',
    'P 63': 'Hexagonal',
    'P 6 2 2': 'Hexagonal',
    'P 61 2 2': 'Hexagonal',
    'P 65 2 2': 'Hexagonal',
    'P 62 2 2': 'Hexagonal',
    'P 64 2 2': 'Hexagonal',
    'P 63 2 2': 'Hexagonal',
    'P 2 3': 'Cubic',
    'F 2 3': 'Cubic',
    'I 2 3': 'Cubic',
    'P 21 3': 'Cubic',
    'I 21 3': 'Cubic',
    'P 4 3 2': 'Cubic',
    'P 42 3 2': 'Cubic',
    'F 4 3 2': 'Cubic',
    'F 41 3 2': 'Cubic',
    'I 4 3 2': 'Cubic',
    'P 43 3 2': 'Cubic',
    'P 41 3 2': 'Cubic',
    'I 41 3 2': 'Cubic'
}

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# ------------------------------------------------
# PDB Query Functions
# ------------------------------------------------
def query_pdb_structures(query):
    try:
        response = requests.post(BASE_SEARCH_URL, json=query, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        logging.error(f"Error while querying PDB: {e}")
        return []

    try:
        data = response.json()
    except ValueError as e:
        logging.error(f"Invalid JSON response: {e}")
        return []

    if "result_set" not in data:
        logging.info("No 'result_set' key in response; returning empty list.")
        return []

    results = []
    for entry in data["result_set"]:
        pdb_id = entry.get("identifier", "")
        score = entry.get("score", "N/A")
        results.append({"pdb_id": pdb_id, "score": score})
    return results

def filter_pdb_results(results, excluded_ids):
    filtered = [r for r in results if r["pdb_id"] not in excluded_ids]
    return filtered

# ------------------------------------------------
# PDB File Space Group Extraction
# ------------------------------------------------
def extract_space_group_from_pdb(pdb_file_path):
    """Extract space group from PDB file header (REMARK 290)"""
    try:
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith('REMARK 290') and 'SPACE GROUP:' in line:
                    space_group = line.split('SPACE GROUP:')[1].strip()
                    return space_group
                elif line.startswith('REMARK') and 'SPACE GROUP' in line and ':' in line:
                    match = re.search(r'SPACE GROUP\s*:\s*(.+)', line)
                    if match:
                        return match.group(1).strip()
    except Exception as e:
        logging.debug(f"Could not extract space group from {pdb_file_path}: {e}")
    return None

# ------------------------------------------------
# Resolution Classification
# ------------------------------------------------
def classify_resolution(resolution):
    """Classify resolution into predefined bins"""
    if resolution is None:
        return None
    
    try:
        res_float = float(resolution)
        for bin_name, test_fn in RESOLUTION_BINS.items():
            if test_fn(res_float):
                return bin_name
    except (ValueError, TypeError):
        pass
    return None

# ------------------------------------------------
# Enhanced B-Factor Extraction Functions
# ------------------------------------------------
def extract_bfactors_from_pdb(pdb_file_path):
    """Extract B-factors from ATOM and HETATM records in PDB file"""
    try:
        b_factors = []
        in_model_1 = True
        
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    model_num = int(line.split()[1]) if len(line.split()) > 1 else 1
                    in_model_1 = (model_num == 1)
                elif line.startswith('ENDMDL'):
                    if in_model_1:
                        break
                elif in_model_1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                    try:
                        altloc = line[16] if len(line) > 16 else ' '
                        if altloc in [' ', 'A']:
                            bfactor_str = line[60:66].strip()
                            if bfactor_str:
                                b_factor = float(bfactor_str)
                                b_factors.append(b_factor)
                    except (ValueError, IndexError):
                        continue
        
        if not b_factors:
            return None
        
        stats = {
            'mean_bfactor': round(statistics.mean(b_factors), 2),
            'min_bfactor': round(min(b_factors), 2),
            'max_bfactor': round(max(b_factors), 2),
            'bfactor_range': round(max(b_factors) - min(b_factors), 2),
            'num_atoms': len(b_factors)
        }
        
        return stats
        
    except Exception as e:
        logging.debug(f"Could not extract B-factors from {pdb_file_path}: {e}")
        return None

def get_bfactor_data(pdb_id, output_dir=OUTPUT_DIR):
    """Get B-factor data for a PDB ID if the file exists"""
    pdb_file = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    if os.path.exists(pdb_file):
        return extract_bfactors_from_pdb(pdb_file)
    return None

# ------------------------------------------------
# VDR Chain Identification Functions
# ------------------------------------------------
def identify_vdr_chains_from_pdb(pdb_file_path):
    """Identify chains containing VDR by parsing COMPND records"""
    try:
        vdr_chains = []
        current_mol_id = None
        current_chains = []
        current_molecule = ""
        
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith('COMPND'):
                    if 'MOL_ID:' in line:
                        if current_mol_id and any(vdr_id.lower() in current_molecule.lower() 
                                                for vdr_id in VDR_IDENTIFIERS):
                            vdr_chains.extend(current_chains)
                        
                        current_mol_id = line.split('MOL_ID:')[1].strip().rstrip(';')
                        current_chains = []
                        current_molecule = ""
                    elif 'MOLECULE:' in line:
                        current_molecule = line.split('MOLECULE:')[1].strip().rstrip(';')
                    elif 'CHAIN:' in line:
                        chains_str = line.split('CHAIN:')[1].strip().rstrip(';')
                        current_chains = [c.strip() for c in chains_str.split(',')]
            
            if current_mol_id and any(vdr_id.lower() in current_molecule.lower() 
                                    for vdr_id in VDR_IDENTIFIERS):
                vdr_chains.extend(current_chains)
        
        return vdr_chains
        
    except Exception as e:
        logging.debug(f"Could not identify VDR chains from {pdb_file_path}: {e}")
        return []

def extract_vdr_organism_from_pdb(pdb_file_path):
    """Extract organism for VDR chains from PDB COMPND records"""
    try:
        vdr_chains = identify_vdr_chains_from_pdb(pdb_file_path)
        if not vdr_chains:
            return None
        
        with open(pdb_file_path, 'r') as f:
            content = f.read()
            
        for species in SPECIES_LIST:
            if species.lower() in content.lower():
                return species
        
        if any(keyword in content.lower() for keyword in ["human", "humans"]):
            return "Homo sapiens"
        elif any(keyword in content.lower() for keyword in ["rat", "rats"]):
            return "Rattus norvegicus"
        elif any(keyword in content.lower() for keyword in ["zebrafish", "zebra fish", "danio"]):
            return "Danio rerio"
        elif any(keyword in content.lower() for keyword in ["lamprey", "sea lamprey", "petromyzon"]):
            return "Petromyzon marinus"
            
    except Exception as e:
        logging.debug(f"Could not extract VDR organism from {pdb_file_path}: {e}")
    
    return None

# ------------------------------------------------
# API-based Crystallization Data Functions
# ------------------------------------------------
def get_crystallization_data_from_api(pdb_id: str, max_retries: int = MAX_RETRIES) -> Dict[str, Any]:
    """Fetch crystallization data from REST API"""
    url = f"{BASE_METADATA_URL}{pdb_id}"
    
    for attempt in range(max_retries + 1):
        try:
            if attempt > 0:
                delay = RETRY_DELAY * (2 ** (attempt - 1))
                logging.info(f"Retrying crystallization data for {pdb_id} after {delay} seconds...")
                time.sleep(delay)
            else:
                time.sleep(REQUEST_DELAY)
            
            resp = requests.get(url, timeout=15)
            
            if resp.status_code == 429:
                retry_after = resp.headers.get('Retry-After')
                if retry_after:
                    wait_time = int(retry_after)
                    logging.warning(f"Rate limited for {pdb_id}. Waiting {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    logging.warning(f"Rate limited for {pdb_id}. Using exponential backoff...")
                continue
            
            resp.raise_for_status()
            
        except requests.exceptions.RequestException as e:
            logging.error(f"Error fetching crystallization data for {pdb_id} (attempt {attempt + 1}): {e}")
            if attempt == max_retries:
                return {}
            continue

        try:
            data = resp.json()
        except ValueError as e:
            logging.error(f"Invalid JSON for {pdb_id}: {e}")
            return {}

        crystal_info = {}
        
        entry_info = data.get("rcsb_entry_info", {})
        if "resolution_combined" in entry_info:
            res_list = entry_info["resolution_combined"]
            if isinstance(res_list, list) and res_list:
                crystal_info['resolution'] = res_list[0]
        else:
            refine = data.get("refine", [])
            if isinstance(refine, list) and refine:
                for refine_entry in refine:
                    if isinstance(refine_entry, dict) and "ls_d_res_high" in refine_entry:
                        crystal_info['resolution'] = refine_entry["ls_d_res_high"]
                        break
        
        symmetry = data.get("symmetry", {})
        if symmetry and isinstance(symmetry, dict):
            space_group = symmetry.get("space_group_name_hm", "").strip()
            if space_group:
                crystal_info['space_group'] = space_group
                crystal_system = CRYSTAL_SYSTEMS.get(space_group, "Unknown")
                crystal_info['crystal_system'] = crystal_system
            else:
                crystal_info['space_group'] = ""
                crystal_info['crystal_system'] = "Unknown"
        else:
            crystal_info['space_group'] = ""
            crystal_info['crystal_system'] = "Unknown"
        
        cell = data.get("cell", {})
        if cell and isinstance(cell, dict):
            crystal_info.update({
                'unit_cell_a': cell.get("length_a"),
                'unit_cell_b': cell.get("length_b"), 
                'unit_cell_c': cell.get("length_c"),
                'unit_cell_alpha': cell.get("angle_alpha"),
                'unit_cell_beta': cell.get("angle_beta"),
                'unit_cell_gamma': cell.get("angle_gamma")
            })
        
        exptl_crystal_grow = data.get("exptl_crystal_grow", [])
        if exptl_crystal_grow and isinstance(exptl_crystal_grow, list):
            grow_info = exptl_crystal_grow[0]
            crystal_info['crystal_method'] = grow_info.get("method", "")
        
        return crystal_info
    
    return {}

# ------------------------------------------------
# VDR Organism Extraction from Web Scraping
# ------------------------------------------------
def extract_vdr_organism_from_table(html_content: str, pdb_id: str) -> Optional[str]:
    """Parse macromolecules table to find VDR organism"""
    try:
        soup = BeautifulSoup(html_content, 'html.parser')
        
        tables = soup.find_all('table')
        
        for table_idx, table in enumerate(tables):
            headers = table.find_all('th')
            header_text = ' '.join([th.get_text().strip().lower() for th in headers])
            
            if any(keyword in header_text for keyword in ['molecule', 'organism', 'chains']):
                rows = table.find_all('tr')
                
                for row_idx, row in enumerate(rows):
                    cells = row.find_all(['td', 'th'])
                    if len(cells) >= 3:
                        row_text = ' '.join([cell.get_text().strip() for cell in cells])
                        row_text_lower = row_text.lower()
                        
                        vdr_identifiers_lower = [vdr_id.lower() for vdr_id in VDR_IDENTIFIERS]
                        if any(vdr_id in row_text_lower for vdr_id in vdr_identifiers_lower):
                            organism = extract_organism_from_row(cells, pdb_id)
                            if organism:
                                return organism
        
        return None
        
    except Exception as e:
        logging.error(f"Error parsing macromolecules table for {pdb_id}: {e}")
        return None

def extract_organism_from_row(cells, pdb_id):
    """Extract organism from table row containing VDR information"""
    cell_texts = [cell.get_text().strip() for cell in cells]
    
    for i, cell_text in enumerate(cell_texts):
        for species in SPECIES_LIST:
            if species.lower() == cell_text.lower().strip():
                return species
        
        for species in SPECIES_LIST:
            if species.lower() in cell_text.lower():
                return species
        
        cell_lower = cell_text.lower().strip()
        if cell_lower in ["human", "humans"] or "homo sapiens" in cell_lower:
            return "Homo sapiens"
        elif cell_lower in ["rat", "rats"] or "rattus norvegicus" in cell_lower:
            return "Rattus norvegicus"
        elif any(name in cell_lower for name in ["zebrafish", "zebra fish", "danio"]):
            return "Danio rerio"
        elif any(name in cell_lower for name in ["lamprey", "sea lamprey", "petromyzon"]):
            return "Petromyzon marinus"
    
    return None

def extract_fallback_organism(html_content, pdb_id):
    """Fallback method to extract organism"""
    metadata_section = html_content[:5000]
    
    for species in SPECIES_LIST:
        if species.lower() in metadata_section.lower():
            return species
    
    if any(keyword in metadata_section.lower() for keyword in ["human", "humans"]):
        return "Homo sapiens"
    elif any(keyword in metadata_section.lower() for keyword in ["rat", "rats"]):
        return "Rattus norvegicus"
    elif any(keyword in metadata_section.lower() for keyword in ["zebrafish", "zebra fish", "danio"]):
        return "Danio rerio"
    elif any(keyword in metadata_section.lower() for keyword in ["lamprey", "sea lamprey", "petromyzon"]):
        return "Petromyzon marinus"
    
    return None

def scrape_structure_info(pdb_id: str, max_retries: int = MAX_RETRIES) -> Tuple[Optional[str], Optional[str]]:
    """Scrape structure page for name and VDR organism"""
    url = f"{BASE_STRUCTURE_URL}{pdb_id}"
    
    for attempt in range(max_retries + 1):
        try:
            if attempt > 0:
                delay = RETRY_DELAY * (2 ** (attempt - 1))
                time.sleep(delay)
            else:
                time.sleep(REQUEST_DELAY)
            
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            
            resp = requests.get(url, timeout=15, headers=headers)
            
            if resp.status_code == 429:
                retry_after = resp.headers.get('Retry-After')
                if retry_after:
                    wait_time = int(retry_after)
                    time.sleep(wait_time)
                else:
                    continue
            
            resp.raise_for_status()
            
        except requests.exceptions.RequestException as e:
            if attempt == max_retries:
                return None, None
            continue

        try:
            html_content = resp.text
            
            structure_name = None
            title_match = re.search(rf'<title[^>]*>.*?{pdb_id}[:\-\s]+([^|<]+)', html_content, re.IGNORECASE | re.DOTALL)
            if title_match:
                structure_name = title_match.group(1).strip()
                structure_name = re.sub(r'\s*\|\s*RCSB.*$', '', structure_name)
                structure_name = re.sub(r'\s*-\s*RCSB.*$', '', structure_name)
                structure_name = structure_name.strip()
            
            vdr_organism = extract_vdr_organism_from_table(html_content, pdb_id)
            
            if vdr_organism:
                return structure_name, vdr_organism
            
            fallback_organism = extract_fallback_organism(html_content, pdb_id)
            
            if fallback_organism:
                return structure_name, f"{fallback_organism} (uncertain)"
            
            return structure_name, None
            
        except Exception as e:
            logging.error(f"Error parsing HTML for {pdb_id}: {e}")
            return None, None
    
    return None, None

# ------------------------------------------------
# Complete Structure Information Assembly
# ------------------------------------------------
def get_complete_structure_info(pdb_id: str) -> Dict[str, Any]:
    """Get complete information for a PDB structure"""
    crystal_data = get_crystallization_data_from_api(pdb_id)
    
    structure_name, organism = scrape_structure_info(pdb_id)
    
    pdb_file_path = os.path.join(OUTPUT_DIR, f"{pdb_id}.pdb")
    if os.path.exists(pdb_file_path):
        pdb_space_group = extract_space_group_from_pdb(pdb_file_path)
        if pdb_space_group and not crystal_data.get('space_group'):
            crystal_data['space_group'] = pdb_space_group
            crystal_data['crystal_system'] = CRYSTAL_SYSTEMS.get(pdb_space_group, "Unknown")
        
        pdb_organism = extract_vdr_organism_from_pdb(pdb_file_path)
        if pdb_organism and not organism:
            organism = pdb_organism
    
    bfactor_data = get_bfactor_data(pdb_id)
    
    structure_info = {
        'pdb_id': pdb_id,
        'structure_name': structure_name,
        'organism': organism,
    }
    
    structure_info.update(crystal_data)
    
    if structure_info.get('resolution'):
        structure_info['resolution_group'] = classify_resolution(structure_info['resolution'])
    
    if bfactor_data:
        structure_info.update(bfactor_data)
    
    return structure_info

# ------------------------------------------------
# Structure Processing and CSV Output
# ------------------------------------------------
def process_and_save_structures(pdb_ids):
    """Process all PDB IDs and save enhanced data to CSV"""
    all_structures = []
    total_ids = len(pdb_ids)
    
    logging.info(f"Processing {total_ids} PDB IDs with crystallization data...")
    logging.info("Extracting: VDR organism, resolution, space group, unit cell, B-factors")
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        logging.info(f"Processing {pdb_id} ({i}/{total_ids})...")
        
        structure_info = get_complete_structure_info(pdb_id)
        all_structures.append(structure_info)
        
        name = structure_info.get('structure_name') or 'Unknown'
        org = structure_info.get('organism') or 'Unknown'
        res = structure_info.get('resolution', 'Unknown')
        res_group = structure_info.get('resolution_group', 'Unknown')
        sg = structure_info.get('space_group', 'Unknown')
        cs = structure_info.get('crystal_system', 'Unknown')
        bfactor = structure_info.get('mean_bfactor', 'N/A')
        
        logging.info(f"  → Name: {name[:60]}{'...' if len(name) > 60 else ''}")
        logging.info(f"  → VDR Organism: {org}")
        logging.info(f"  → Resolution: {res} Å ({res_group})")
        logging.info(f"  → Space Group: {sg} ({cs})")
        if bfactor != 'N/A':
            logging.info(f"  → Mean B-factor: {bfactor} Ų")
    
    save_to_csv(all_structures)
    analyze_space_groups()
    create_vdr_structure_table()
    
    return all_structures

def save_to_csv(structures):
    """Save structure information to CSV"""
    fieldnames = [
        'pdb_id', 'structure_name', 'organism', 'resolution', 'resolution_group',
        'space_group', 'crystal_system', 
        'unit_cell_a', 'unit_cell_b', 'unit_cell_c',
        'unit_cell_alpha', 'unit_cell_beta', 'unit_cell_gamma',
        'crystal_method', 'mean_bfactor', 'min_bfactor', 'max_bfactor', 'bfactor_range'
    ]
    
    with open(CSV_FILENAME, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        for structure in structures:
            row_data = {}
            for field in fieldnames:
                value = structure.get(field)
                if value is not None:
                    row_data[field] = value
                else:
                    row_data[field] = ''
            writer.writerow(row_data)
    
    logging.info(f"Data saved to {CSV_FILENAME}")
    print(f"\n✅ Data saved to {CSV_FILENAME}")

# ------------------------------------------------
# Analysis Functions
# ------------------------------------------------
def analyze_space_groups():
    """Analyze space group distribution with B-factor analysis"""
    if not os.path.exists(CSV_FILENAME):
        print(f"❌ Error: {CSV_FILENAME} not found!")
        return
    
    print(f"\n🔍 SPACE GROUP ANALYSIS")
    print("="*50)
    
    try:
        df = pd.read_csv(CSV_FILENAME)
        
        df_with_sg = df[df['space_group'].notna() & (df['space_group'] != '')]
        
        print(f"📊 Total structures: {len(df)}")
        print(f"📊 Structures with space group data: {len(df_with_sg)}")
        
        if len(df_with_sg) == 0:
            print("❌ No space group data available for analysis")
            return
        
        sg_counts = df_with_sg['space_group'].value_counts()
        print(f"\n🏗️  SPACE GROUP DISTRIBUTION:")
        print("-" * 30)
        for sg, count in sg_counts.head(10).items():
            percentage = (count / len(df_with_sg)) * 100
            crystal_sys = CRYSTAL_SYSTEMS.get(sg, "Unknown")
            print(f"  {sg:<15} | {crystal_sys:<12} | {count:>3} ({percentage:>5.1f}%)")
        
        df_with_sg['crystal_system'] = df_with_sg['space_group'].map(
            lambda x: CRYSTAL_SYSTEMS.get(x, "Unknown")
        )
        cs_counts = df_with_sg['crystal_system'].value_counts()
        
        print(f"\n🔬 CRYSTAL SYSTEM DISTRIBUTION:")
        print("-" * 30)
        for cs, count in cs_counts.items():
            percentage = (count / len(df_with_sg)) * 100
            print(f"  {cs:<15} | {count:>3} ({percentage:>5.1f}%)")
        
        print(f"\n🧬 SPACE GROUPS BY VDR ORGANISM:")
        print("-" * 40)
        
        organism_sg = df_with_sg.groupby(['organism', 'space_group']).size().unstack(fill_value=0)
        for organism in organism_sg.index:
            if organism and organism != 'Unknown' and not organism.startswith('Unknown'):
                print(f"\n  {organism}:")
                org_sgs = organism_sg.loc[organism]
                org_sgs_nonzero = org_sgs[org_sgs > 0].sort_values(ascending=False)
                for sg, count in org_sgs_nonzero.items():
                    crystal_sys = CRYSTAL_SYSTEMS.get(sg, "Unknown")
                    print(f"    {sg:<15} ({crystal_sys:<12}): {count}")
        
        if 'resolution' in df_with_sg.columns:
            analyze_space_group_vs_resolution(df_with_sg)
        
        if 'mean_bfactor' in df_with_sg.columns:
            analyze_space_group_vs_bfactor(df_with_sg)
        
    except Exception as e:
        print(f"❌ Error analyzing space groups: {e}")

def analyze_space_group_vs_resolution(df):
    """Analyze relationship between space groups and resolution"""
    print(f"\n📐 SPACE GROUP vs RESOLUTION ANALYSIS:")
    print("-" * 45)
    
    df_complete = df[(df['space_group'].notna()) & 
                     (df['resolution'].notna()) & 
                     (df['space_group'] != '') & 
                     (df['resolution'] != '')]
    
    if len(df_complete) == 0:
        print("❌ No structures with both space group and resolution data")
        return
    
    df_complete['resolution'] = pd.to_numeric(df_complete['resolution'], errors='coerce')
    df_complete = df_complete[df_complete['resolution'].notna()]
    
    sg_resolution = df_complete.groupby('space_group')['resolution'].agg(['mean', 'std', 'count'])
    sg_resolution = sg_resolution[sg_resolution['count'] >= 2]
    sg_resolution = sg_resolution.sort_values('mean')
    
    print(f"Average resolution by space group (≥2 structures):")
    print(f"{'Space Group':<15} | {'Avg Res (Å)':<12} | {'Std Dev':<8} | {'Count':<5}")
    print("-" * 50)
    
    for sg, data in sg_resolution.iterrows():
        print(f"{sg:<15} | {data['mean']:>8.2f}     | {data['std']:>6.2f} | {data['count']:>3.0f}")

def analyze_space_group_vs_bfactor(df):
    """Analyze relationship between space groups and B-factors"""
    print(f"\n🌡️  SPACE GROUP vs B-FACTOR ANALYSIS:")
    print("-" * 45)
    
    df_complete = df[(df['space_group'].notna()) & 
                     (df['mean_bfactor'].notna()) & 
                     (df['space_group'] != '') & 
                     (df['mean_bfactor'] != '')]
    
    if len(df_complete) == 0:
        print("❌ No structures with both space group and B-factor data")
        print("💡 Download PDB files to enable B-factor analysis")
        return
    
    df_complete['mean_bfactor'] = pd.to_numeric(df_complete['mean_bfactor'], errors='coerce')
    df_complete = df_complete[df_complete['mean_bfactor'].notna()]
    
    print(f"📊 Analyzing {len(df_complete)} structures with B-factor data")
    
    sg_bfactor = df_complete.groupby('space_group')['mean_bfactor'].agg(['mean', 'std', 'count'])
    sg_bfactor = sg_bfactor[sg_bfactor['count'] >= 2]
    sg_bfactor = sg_bfactor.sort_values('mean')
    
    print(f"\nMean B-factor by space group (≥2 structures):")
    print(f"{'Space Group':<15} | {'Mean B-fact':<11} | {'Std Dev':<8} | {'Count':<5}")
    print("-" * 50)
    
    for sg, data in sg_bfactor.iterrows():
        crystal_sys = CRYSTAL_SYSTEMS.get(sg, "Unknown")
        print(f"{sg:<15} | {data['mean']:>8.1f} Ų  | {data['std']:>6.1f} | {data['count']:>3.0f}")
    
    df_complete['crystal_system'] = df_complete['space_group'].map(
        lambda x: CRYSTAL_SYSTEMS.get(x, "Unknown")
    )
    
    cs_bfactor = df_complete.groupby('crystal_system')['mean_bfactor'].agg(['mean', 'std', 'count'])
    cs_bfactor = cs_bfactor.sort_values('mean')
    
    print(f"\nMean B-factor by crystal system:")
    print(f"{'Crystal System':<15} | {'Mean B-fact':<11} | {'Std Dev':<8} | {'Count':<5}")
    print("-" * 50)
    
    for cs, data in cs_bfactor.iterrows():
        print(f"{cs:<15} | {data['mean']:>8.1f} Ų  | {data['std']:>6.1f} | {data['count']:>3.0f}")
    
    print(f"\n💡 INTERPRETATION:")
    print(f"Higher B-factors = More atomic flexibility/disorder")
    print(f"Lower B-factors = More rigid/ordered structure")
    print(f"This shows how different crystal packings affect VDR flexibility!")

def create_vdr_structure_table():
    """Create VDR structure distribution table by species and resolution"""
    print(f"\n📊 Creating VDR structure distribution table...")
    
    species_display_names = {
        "Homo sapiens": "Human",
        "Rattus norvegicus": "Rats", 
        "Danio rerio": "Zebra Fish",
        "Petromyzon marinus": "Sea Lamprey"
    }
    
    table_data = defaultdict(lambda: defaultdict(int))
    
    if not os.path.exists(CSV_FILENAME):
        print(f"❌ Error: {CSV_FILENAME} not found!")
        return
    
    try:
        with open(CSV_FILENAME, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            total_structures = 0
            processed_structures = 0
            
            for row in reader:
                total_structures += 1
                pdb_id = row['pdb_id']
                organism = row['organism'].strip()
                resolution_str = row['resolution'].strip()
                
                if not resolution_str or resolution_str == '':
                    continue
                
                try:
                    resolution = float(resolution_str)
                except ValueError:
                    continue
                
                clean_organism = organism.replace(" (uncertain)", "")
                
                resolution_bin = None
                for bin_name, test_fn in RESOLUTION_BINS.items():
                    if test_fn(resolution):
                        resolution_bin = bin_name
                        break
                
                if resolution_bin is None:
                    continue
                
                if clean_organism in species_display_names:
                    display_name = species_display_names[clean_organism]
                    table_data[resolution_bin][display_name] += 1
                    processed_structures += 1
            
            print(f"📊 Processed {processed_structures} out of {total_structures} structures from {CSV_FILENAME}")
    
    except Exception as e:
        print(f"❌ Error reading {CSV_FILENAME}: {e}")
        return
    
    print_vdr_table(table_data, species_display_names)

def print_vdr_table(table_data, species_display_names):
    """Print the VDR structure table"""
    species_order = ["Human", "Rats", "Zebra Fish", "Sea Lamprey"]
    resolution_order = ["1-2 Å", "2-2.5 Å", "2.5+ Å"]
    
    max_species_width = max(len(species) for species in species_order)
    col_width = max(max_species_width, 8)
    
    print(f"\n{'':>12}", end="")
    for species in species_order:
        print(f"{species:>{col_width}}", end="")
    print()
    
    print("-" * (12 + col_width * len(species_order)))
    
    totals = defaultdict(int)
    for res_bin in resolution_order:
        print(f"{res_bin:>12}", end="")
        for species in species_order:
            count = table_data[res_bin][species]
            print(f"{count if count > 0 else '-':>{col_width}}", end="")
            totals[species] += count
        print()
    
    print("-" * (12 + col_width * len(species_order)))
    
    print(f"{'Total:':>12}", end="")
    grand_total = 0
    for species in species_order:
        total = totals[species]
        print(f"{total if total > 0 else '-':>{col_width}}", end="")
        grand_total += total
    print()
    
    print(f"\n📈 Grand Total: {grand_total} VDR structures")
    
    print(f"\n📊 VDR Statistics:")
    for species in species_order:
        if totals[species] > 0:
            percentage = (totals[species] / grand_total) * 100
            print(f"  VDR from {species}: {totals[species]} structures ({percentage:.1f}%)")

# ------------------------------------------------
# File Download Functions
# ------------------------------------------------
def download_structures(pdb_ids, output_dir=OUTPUT_DIR):
    """Download structure files for given PDB IDs"""
    os.makedirs(output_dir, exist_ok=True)
    downloaded = []

    print(f"\n📥 Downloading {len(pdb_ids)} structure files...")
    print("💡 PDB files are needed for B-factor analysis")

    for i, pdb_id in enumerate(pdb_ids, 1):
        pdb_id = pdb_id.upper()
        
        if i % 20 == 0:
            print(f"   Downloaded {i}/{len(pdb_ids)} files...")

        pdb_url = f"{BASE_DOWNLOAD_URL}{pdb_id}.pdb"
        try:
            response = requests.get(pdb_url, timeout=10)
        except requests.exceptions.RequestException as e:
            logging.error(f"Network error when downloading {pdb_id}.pdb: {e}")
            response = None

        if response is not None and response.status_code == 200:
            file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
            try:
                with open(file_path, "w") as fout:
                    fout.write(response.text)
                downloaded.append((pdb_id, ".pdb"))
                continue
            except IOError as e:
                logging.error(f"Error saving {pdb_id}.pdb: {e}")
        else:
            cif_url = f"{BASE_DOWNLOAD_URL}{pdb_id}.cif"
            try:
                cif_response = requests.get(cif_url, timeout=10)
            except requests.exceptions.RequestException as e:
                continue

            if cif_response is not None and cif_response.status_code == 200:
                file_path = os.path.join(output_dir, f"{pdb_id}.cif")
                try:
                    with open(file_path, "w") as fout:
                        fout.write(cif_response.text)
                    downloaded.append((pdb_id, ".cif"))
                    continue
                except IOError as e:
                    continue

    return downloaded

# ------------------------------------------------
# Main Execution
# ------------------------------------------------
def main():
    """Enhanced main function with space group and B-factor analysis"""
    print("🔍 Enhanced VDR Structure Extraction Tool")
    print("="*70)
    print("Extracts: VDR organisms, resolution, space groups, unit cells, B-factors")
    print("Focus: Crystallization effects on VDR protein conformation\n")
    
    query = {
        "query": {
            "type": "group",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": "VDR"},
                },
                {
                    "type": "group",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                                "operator": "in",
                                "value": SPECIES_LIST,
                            },
                        }
                    ],
                    "logical_operator": "or",
                    "label": "rcsb_entity_source_organism.ncbi_scientific_name",
                },
            ],
            "logical_operator": "and",
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 200},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined",
        },
    }
    
    logging.info("Querying PDB database...")
    all_results = query_pdb_structures(query)
    if not all_results:
        logging.info("No results returned from PDB query. Exiting.")
        return

    filtered = filter_pdb_results(all_results, NOT_ALLOWED_PDB_IDS)
    if not filtered:
        logging.info("No structures found after filtering. Exiting.")
        return

    pdb_ids_to_process = [entry["pdb_id"] for entry in filtered]
    print(f"\n📊 Processing {len(pdb_ids_to_process)} structures...")
    
    download_choice = input(f"\n❓ Download structure files for B-factor analysis? (y/N): ").strip().lower()
    if download_choice in ['y', 'yes']:
        downloaded_info = download_structures(pdb_ids_to_process, OUTPUT_DIR)
        pdb_count = len([x for x in downloaded_info if x[1] == '.pdb'])
        print(f"\n✅ Downloaded {len(downloaded_info)} files ({pdb_count} PDB files for B-factor analysis)")
    
    all_structures = process_and_save_structures(pdb_ids_to_process)
    
    print(f"\n🎉 Processing complete!")
    print(f"📋 Check '{CSV_FILENAME}' for detailed crystallization and B-factor data")
    print(f"🔬 Analysis shows how crystallization affects VDR conformation and flexibility")

if __name__ == "__main__":
    main()