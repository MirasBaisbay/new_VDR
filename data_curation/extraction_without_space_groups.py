import os
import logging
import requests
import time
import csv
import re
from typing import Tuple, Optional, Dict, Any
from bs4 import BeautifulSoup
from collections import defaultdict

# -----------------------
# Constants / Configuration
# -----------------------
BASE_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json="
BASE_DOWNLOAD_URL = "https://files.rcsb.org/download/"
BASE_METADATA_URL = "https://data.rcsb.org/rest/v1/core/entry/"
BASE_STRUCTURE_URL = "https://www.rcsb.org/structure/"

# These PDB IDs should be excluded because they are DNA‐binding
NOT_ALLOWED_PDB_IDS = ["1KB2", "1KB4", "1KB6", "1YNW"]

# Where to save downloaded files
OUTPUT_DIR = "pdb_files"
CSV_FILENAME = "vdr_structures.csv"

# Species of interest (in priority order for display)
SPECIES_LIST = [
    "Homo sapiens",
    "Rattus norvegicus", 
    "Danio rerio",
    "Petromyzon marinus",
]

# Resolution bins (in Å)
RESOLUTION_BINS = {
    "1-2 Å": lambda r: 1.0 <= r < 2.0,
    "2-2.5 Å": lambda r: 2.0 <= r < 2.5,
    "2.5+ Å": lambda r: r >= 2.5
}

# Rate limiting settings
REQUEST_DELAY = 0.5  # Delay between requests in seconds (2 requests per second)
MAX_RETRIES = 3      # Maximum number of retries for failed requests
RETRY_DELAY = 5      # Initial delay for exponential backoff in seconds

# VDR-specific search patterns
VDR_IDENTIFIERS = [
    "vitamin d receptor",
    "vdr",
    "vitamin d3 receptor",
    "calcitriol receptor",
    "vitamin nuclear receptor",
    "vitamin d nuclear receptor"
]

# -------------------
# Simple Logging Setup
# -------------------
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


# ------------------------
# Function: Query PDB API
# ------------------------
def query_pdb_structures(query):
    """
    Send a POST request to the RCSB PDB search service with the given JSON query.
    Returns a list of dicts, each containing 'pdb_id' and 'score'.
    """
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


# --------------------------
# Function: Filter PDB Results
# --------------------------
def filter_pdb_results(results, excluded_ids):
    """
    Remove any entries from `results` whose 'pdb_id' is in `excluded_ids`.
    """
    filtered = [r for r in results if r["pdb_id"] not in excluded_ids]
    return filtered


# ------------------------------------------------
# Function: Get resolution from REST API
# ------------------------------------------------
def get_resolution_from_api(pdb_id: str, max_retries: int = MAX_RETRIES) -> Optional[float]:
    """
    Fetch resolution for a PDB ID via the RCSB REST API.
    Returns resolution value or None on failure.
    """
    url = f"{BASE_METADATA_URL}{pdb_id}"
    
    for attempt in range(max_retries + 1):
        try:
            # Add delay between requests
            if attempt > 0:
                delay = RETRY_DELAY * (2 ** (attempt - 1))
                logging.info(f"Retrying resolution for {pdb_id} after {delay} seconds...")
                time.sleep(delay)
            else:
                time.sleep(REQUEST_DELAY)
            
            resp = requests.get(url, timeout=15)
            
            # Handle rate limiting
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
            logging.error(f"Error fetching resolution for {pdb_id} (attempt {attempt + 1}): {e}")
            if attempt == max_retries:
                return None
            continue

        # Parse JSON response
        try:
            data = resp.json()
        except ValueError as e:
            logging.error(f"Invalid JSON for {pdb_id}: {e}")
            return None

        # Extract resolution from resolution_combined field
        entry_info = data.get("rcsb_entry_info", {})
        if "resolution_combined" in entry_info:
            res_list = entry_info["resolution_combined"]
            if isinstance(res_list, list) and res_list:
                return res_list[0]
        
        # Fallback to other fields if resolution_combined not found
        refine = data.get("refine", [])
        if isinstance(refine, list) and refine:
            for refine_entry in refine:
                if isinstance(refine_entry, dict) and "ls_d_res_high" in refine_entry:
                    return refine_entry["ls_d_res_high"]
        
        return None
    
    return None


# ------------------------------------------------
# Function: Extract VDR organism from macromolecules table
# ------------------------------------------------
def extract_vdr_organism_from_table(html_content: str, pdb_id: str) -> Optional[str]:
    """
    Parse the macromolecules table to find the organism of the VDR protein specifically.
    Returns the organism of the VDR protein, or None if not found.
    """
    try:
        soup = BeautifulSoup(html_content, 'html.parser')
        
        # Method 1: Look for the macromolecules table
        # The table typically has headers like: Molecule | Chains | Sequence Length | Organism | Details | Image
        
        # Find tables that might contain macromolecule information
        tables = soup.find_all('table')
        
        for table_idx, table in enumerate(tables):
            # Check if this table contains macromolecule information
            headers = table.find_all('th')
            header_text = ' '.join([th.get_text().strip().lower() for th in headers])
            
            # Look for tables with relevant headers
            if any(keyword in header_text for keyword in ['molecule', 'organism', 'chains']):
                logging.debug(f"{pdb_id}: Found potential macromolecules table {table_idx} with headers: {header_text}")
                
                # Parse table rows
                rows = table.find_all('tr')
                logging.debug(f"{pdb_id}: Table {table_idx} has {len(rows)} rows")
                
                for row_idx, row in enumerate(rows):
                    cells = row.find_all(['td', 'th'])
                    if len(cells) >= 3:  # Relax requirement - need at least molecule, chains, organism
                        row_text = ' '.join([cell.get_text().strip() for cell in cells])
                        row_text_lower = row_text.lower()
                        
                        logging.debug(f"{pdb_id}: Row {row_idx} ({len(cells)} cells): {row_text[:150]}...")
                        
                        # Check if this row contains VDR - use lowercase for comparison
                        vdr_identifiers_lower = [vdr_id.lower() for vdr_id in VDR_IDENTIFIERS]
                        if any(vdr_id in row_text_lower for vdr_id in vdr_identifiers_lower):
                            logging.debug(f"{pdb_id}: Found VDR row {row_idx}: {row_text}")
                            
                            # Try to extract organism from this row
                            organism = extract_organism_from_row(cells, pdb_id)
                            if organism:
                                logging.info(f"{pdb_id}: Found VDR organism from table: {organism}")
                                return organism
                            else:
                                logging.debug(f"{pdb_id}: VDR row found but no organism extracted")
        
        # Method 2: Alternative table structure - look for entity tables
        # Some pages have different table structures
        entity_divs = soup.find_all('div', {'id': re.compile(r'entity.*', re.IGNORECASE)})
        for div in entity_divs:
            div_text = div.get_text()
            div_text_lower = div_text.lower()
            
            if any(vdr_id.lower() in div_text_lower for vdr_id in VDR_IDENTIFIERS):
                logging.debug(f"{pdb_id}: Found VDR in entity div: {div_text[:100]}...")
                organism = extract_organism_from_text_simple(div_text)
                if organism:
                    logging.info(f"{pdb_id}: Found VDR organism from entity div: {organism}")
                    return organism
        
        # Method 3: Look for structured data in the page
        # Find sections that mention VDR
        for vdr_id in VDR_IDENTIFIERS:
            # Case insensitive search for VDR mentions
            pattern = re.compile(re.escape(vdr_id), re.IGNORECASE)
            matches = pattern.finditer(str(soup))
            
            for match in matches:
                # Get context around the match
                start = max(0, match.start() - 500)
                end = min(len(str(soup)), match.end() + 500)
                context = str(soup)[start:end]
                
                # Look for organism information in this context
                organism = extract_organism_from_context(context, pdb_id)
                if organism:
                    logging.info(f"{pdb_id}: Found VDR organism from context: {organism}")
                    return organism
        
        logging.warning(f"{pdb_id}: Could not find VDR organism in macromolecules table")
        return None
        
    except Exception as e:
        logging.error(f"Error parsing macromolecules table for {pdb_id}: {e}")
        return None


def extract_organism_from_row(cells, pdb_id):
    """
    Extract organism from a table row that contains VDR information.
    """
    # Convert cells to text
    cell_texts = [cell.get_text().strip() for cell in cells]
    
    logging.debug(f"{pdb_id}: Analyzing {len(cell_texts)} cells: {cell_texts}")
    
    # Look for organism names in the cells
    for i, cell_text in enumerate(cell_texts):
        # Check if this cell contains a known species (exact match)
        for species in SPECIES_LIST:
            if species.lower() == cell_text.lower().strip():
                logging.debug(f"{pdb_id}: Found exact species match {species} in cell {i}: '{cell_text}'")
                return species
        
        # Check if species name is contained within the cell text
        for species in SPECIES_LIST:
            if species.lower() in cell_text.lower():
                logging.debug(f"{pdb_id}: Found species {species} contained in cell {i}: '{cell_text}'")
                return species
        
        # Also check for common names
        cell_lower = cell_text.lower().strip()
        if cell_lower in ["human", "humans"] or "homo sapiens" in cell_lower:
            logging.debug(f"{pdb_id}: Found human indicator in cell {i}: '{cell_text}'")
            return "Homo sapiens"
        elif cell_lower in ["rat", "rats"] or "rattus norvegicus" in cell_lower:
            logging.debug(f"{pdb_id}: Found rat indicator in cell {i}: '{cell_text}'")
            return "Rattus norvegicus"
        elif any(name in cell_lower for name in ["zebrafish", "zebra fish", "danio"]):
            logging.debug(f"{pdb_id}: Found zebrafish indicator in cell {i}: '{cell_text}'")
            return "Danio rerio"
        elif any(name in cell_lower for name in ["lamprey", "sea lamprey", "petromyzon"]):
            logging.debug(f"{pdb_id}: Found lamprey indicator in cell {i}: '{cell_text}'")
            return "Petromyzon marinus"
    
    logging.debug(f"{pdb_id}: No organism found in any cell")
    return None


def extract_organism_from_text_simple(text):
    """
    Simple organism extraction from text.
    """
    text_lower = text.lower()
    
    # Check for exact species matches first
    for species in SPECIES_LIST:
        if species.lower() in text_lower:
            return species
    
    # Check for common names
    if any(keyword in text_lower for keyword in ["human", "humans"]):
        return "Homo sapiens"
    elif any(keyword in text_lower for keyword in ["rat", "rats"]):
        return "Rattus norvegicus"
    elif any(keyword in text_lower for keyword in ["zebrafish", "zebra fish", "danio"]):
        return "Danio rerio"
    elif any(keyword in text_lower for keyword in ["lamprey", "sea lamprey", "petromyzon"]):
        return "Petromyzon marinus"
    
    return None


def extract_organism_from_context(context, pdb_id):
    """
    Extract organism from HTML context around VDR mentions.
    """
    # Remove HTML tags for easier parsing
    clean_text = re.sub(r'<[^>]+>', ' ', context)
    clean_text = re.sub(r'\s+', ' ', clean_text).strip()
    
    # Look for organism patterns near VDR
    for species in SPECIES_LIST:
        if species.lower() in clean_text.lower():
            # Check if this species appears close to VDR mention
            vdr_pos = -1
            for vdr_id in VDR_IDENTIFIERS:
                pos = clean_text.lower().find(vdr_id)
                if pos != -1:
                    vdr_pos = pos
                    break
            
            if vdr_pos != -1:
                species_pos = clean_text.lower().find(species.lower())
                if species_pos != -1 and abs(species_pos - vdr_pos) < 200:  # Within 200 characters
                    logging.debug(f"{pdb_id}: Found {species} near VDR in context")
                    return species
    
    return None


# ------------------------------------------------
# Function: Scrape structure page for name and VDR organism
# ------------------------------------------------
def scrape_structure_info(pdb_id: str, max_retries: int = MAX_RETRIES) -> Tuple[Optional[str], Optional[str]]:
    """
    Scrape the RCSB structure page to get structure name and VDR-specific organism.
    Focus on finding the organism of the VDR protein specifically.
    """
    url = f"{BASE_STRUCTURE_URL}{pdb_id}"
    
    for attempt in range(max_retries + 1):
        try:
            if attempt > 0:
                delay = RETRY_DELAY * (2 ** (attempt - 1))
                logging.info(f"Retrying scrape for {pdb_id} after {delay} seconds...")
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
                    logging.warning(f"Rate limited for {pdb_id}. Waiting {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    logging.warning(f"Rate limited for {pdb_id}. Using exponential backoff...")
                continue
            
            resp.raise_for_status()
            
        except requests.exceptions.RequestException as e:
            logging.error(f"Error scraping {pdb_id} (attempt {attempt + 1}): {e}")
            if attempt == max_retries:
                return None, None
            continue

        try:
            html_content = resp.text
            
            # Extract structure name from title
            structure_name = None
            title_match = re.search(rf'<title[^>]*>.*?{pdb_id}[:\-\s]+([^|<]+)', html_content, re.IGNORECASE | re.DOTALL)
            if title_match:
                structure_name = title_match.group(1).strip()
                # Clean up common title artifacts
                structure_name = re.sub(r'\s*\|\s*RCSB.*$', '', structure_name)
                structure_name = re.sub(r'\s*-\s*RCSB.*$', '', structure_name)
                structure_name = structure_name.strip()
            
            # NEW APPROACH: Extract VDR-specific organism from macromolecules table
            vdr_organism = extract_vdr_organism_from_table(html_content, pdb_id)
            
            if vdr_organism:
                logging.info(f"{pdb_id}: Successfully found VDR organism: {vdr_organism}")
                return structure_name, vdr_organism
            
            # Fallback: If we can't find VDR-specific organism, try the old method
            # but mark it as uncertain
            logging.warning(f"{pdb_id}: Could not find VDR-specific organism, using fallback method")
            fallback_organism = extract_fallback_organism(html_content, pdb_id)
            
            if fallback_organism:
                logging.info(f"{pdb_id}: Using fallback organism: {fallback_organism}")
                return structure_name, f"{fallback_organism} (uncertain)"
            
            logging.warning(f"{pdb_id}: No organism found with any method")
            return structure_name, None
            
        except Exception as e:
            logging.error(f"Error parsing HTML for {pdb_id}: {e}")
            return None, None
    
    return None, None


def extract_fallback_organism(html_content, pdb_id):
    """
    Fallback method to extract organism when VDR-specific extraction fails.
    Uses the old organism extraction logic as a last resort.
    """
    # Look for organisms in the first part of the page (metadata section)
    metadata_section = html_content[:5000]
    
    for species in SPECIES_LIST:
        if species.lower() in metadata_section.lower():
            return species
    
    # Check for common names
    if any(keyword in metadata_section.lower() for keyword in ["human", "humans"]):
        return "Homo sapiens"
    elif any(keyword in metadata_section.lower() for keyword in ["rat", "rats"]):
        return "Rattus norvegicus"
    elif any(keyword in metadata_section.lower() for keyword in ["zebrafish", "zebra fish", "danio"]):
        return "Danio rerio"
    elif any(keyword in metadata_section.lower() for keyword in ["lamprey", "sea lamprey", "petromyzon"]):
        return "Petromyzon marinus"
    
    return None


# ------------------------------------------------
# Function: Get complete structure metadata
# ------------------------------------------------
def get_complete_structure_info(pdb_id: str) -> Dict[str, Any]:
    """
    Get complete information for a PDB structure including name, VDR organism, and resolution.
    """
    logging.info(f"Getting complete info for {pdb_id}...")
    
    # Get resolution from API
    resolution = get_resolution_from_api(pdb_id)
    
    # Get name and VDR-specific organism from web scraping
    structure_name, organism = scrape_structure_info(pdb_id)
    
    return {
        'pdb_id': pdb_id,
        'structure_name': structure_name,
        'organism': organism,
        'resolution': resolution
    }


# ------------------------------------------------------
# Function: Process and save structures
# ------------------------------------------------------
def process_and_save_structures(pdb_ids):
    """
    Process all PDB IDs and save to CSV.
    """
    # Collect all structure information
    all_structures = []
    total_ids = len(pdb_ids)
    
    logging.info(f"Processing {total_ids} PDB IDs...")
    logging.info("Using VDR-specific organism detection")
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        logging.info(f"Processing {pdb_id} ({i}/{total_ids})...")
        
        structure_info = get_complete_structure_info(pdb_id)
        all_structures.append(structure_info)
        
        # Log the results
        name = structure_info['structure_name'] or 'Unknown'
        org = structure_info['organism'] or 'Unknown'
        res = structure_info['resolution'] if structure_info['resolution'] is not None else 'Unknown'
        logging.info(f"  → Name: {name[:60]}{'...' if len(name) > 60 else ''}")
        logging.info(f"  → VDR Organism: {org}")
        logging.info(f"  → Resolution: {res} Å")
    
    # Save to CSV
    save_to_csv(all_structures)
    
    # Categorize by species and resolution (for display)
    categorize_and_display(all_structures)
    
    return all_structures


def save_to_csv(structures):
    """
    Save structure information to CSV file.
    """
    csv_path = CSV_FILENAME
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['pdb_id', 'structure_name', 'organism', 'resolution']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header
        writer.writeheader()
        
        # Write data
        for structure in structures:
            writer.writerow({
                'pdb_id': structure['pdb_id'],
                'structure_name': structure['structure_name'] or '',
                'organism': structure['organism'] or '',
                'resolution': structure['resolution'] if structure['resolution'] is not None else ''
            })
    
    logging.info(f"Data saved to {csv_path}")
    print(f"\n✅ Data saved to {csv_path}")


def categorize_and_display(structures):
    """
    Categorize structures by VDR species and resolution for display.
    """
    # Initialize nested dict: species → resolution bin → list
    categories = {
        species: {bin_label: [] for bin_label in RESOLUTION_BINS.keys()}
        for species in SPECIES_LIST
    }
    
    # Add categories for uncertain and unknown organisms
    categories["Uncertain"] = {bin_label: [] for bin_label in RESOLUTION_BINS.keys()}
    categories["Unknown"] = {bin_label: [] for bin_label in RESOLUTION_BINS.keys()}
    
    for structure in structures:
        org = structure['organism'] or 'Unknown'
        res = structure['resolution']
        pdb_id = structure['pdb_id']
        
        # Determine category
        if "(uncertain)" in org:
            category_org = "Uncertain"
        elif org == 'Unknown':
            category_org = "Unknown"
        elif org in SPECIES_LIST:
            category_org = org
        else:
            category_org = "Unknown"
        
        if res is not None:
            # Determine which resolution bin it falls into
            placed = False
            for bin_label, test_fn in RESOLUTION_BINS.items():
                if test_fn(res):
                    categories[category_org][bin_label].append(pdb_id)
                    placed = True
                    break
            
            if not placed:
                logging.warning(f"{pdb_id}: resolution {res} Å did not fit any bin")
        else:
            logging.warning(f"{pdb_id}: no resolution data available")
    
    # Display categorized results
    print("\n" + "="*50)
    print("CATEGORIZED RESULTS (VDR-specific organisms)")
    print("="*50)
    print("Note: Organisms are identified based on the VDR protein specifically")
    print("'Uncertain' indicates fallback method was used\n")
    
    total_structures = 0
    for species, bins in categories.items():
        species_count = sum(len(id_list) for id_list in bins.values())
        if species_count > 0:
            total_structures += species_count
            print(f"\nVDR Species: {species} ({species_count} structures)")
            for bin_label, id_list in bins.items():
                if id_list:
                    print(f"  {bin_label}: {len(id_list)} structure(s)")
                    # Show first few IDs to avoid cluttering
                    if len(id_list) <= 10:
                        print(f"    IDs: {', '.join(id_list)}")
                    else:
                        print(f"    IDs: {', '.join(id_list[:10])} ... and {len(id_list)-10} more")
    
    print(f"\nTotal categorized structures: {total_structures}")
    
    # Show summary statistics
    print(f"\n" + "="*30)
    print("SUMMARY STATISTICS")
    print("="*30)
    for species in SPECIES_LIST:
        if species in categories:
            count = sum(len(id_list) for id_list in categories[species].values())
            if count > 0:
                print(f"VDR from {species}: {count} structures")
    
    uncertain_count = sum(len(id_list) for id_list in categories.get("Uncertain", {}).values())
    if uncertain_count > 0:
        print(f"Uncertain VDR species: {uncertain_count} structures")
    
    unknown_count = sum(len(id_list) for id_list in categories.get("Unknown", {}).values())
    if unknown_count > 0:
        print(f"Unknown VDR species: {unknown_count} structures")


# ------------------------------------------------------
# Function: Download Structures
# ------------------------------------------------------
def download_structures(pdb_ids, output_dir=OUTPUT_DIR):
    """
    Download structure files for the given PDB IDs.
    """
    os.makedirs(output_dir, exist_ok=True)
    downloaded = []

    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()

        # First try .pdb
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
                logging.info(f"Downloaded {pdb_id}.pdb → {file_path}")
                downloaded.append((pdb_id, ".pdb"))
                continue
            except IOError as e:
                logging.error(f"Error saving {pdb_id}.pdb: {e}")
        else:
            status = getattr(response, "status_code", "N/A")
            logging.info(f"{pdb_id}.pdb not available (status: {status}). Trying .cif...")

        # Fallback to .cif
        cif_url = f"{BASE_DOWNLOAD_URL}{pdb_id}.cif"
        try:
            cif_response = requests.get(cif_url, timeout=10)
        except requests.exceptions.RequestException as e:
            logging.error(f"Network error when downloading {pdb_id}.cif: {e}")
            cif_response = None

        if cif_response is not None and cif_response.status_code == 200:
            file_path = os.path.join(output_dir, f"{pdb_id}.cif")
            try:
                with open(file_path, "w") as fout:
                    fout.write(cif_response.text)
                logging.info(f"Downloaded {pdb_id}.cif → {file_path}")
                downloaded.append((pdb_id, ".cif"))
                continue
            except IOError as e:
                logging.error(f"Error saving {pdb_id}.cif: {e}")
        else:
            cif_status = getattr(cif_response, "status_code", "N/A")
            logging.error(
                f"Failed to download both {pdb_id}.pdb (status: {status}) "
                f"and {pdb_id}.cif (status: {cif_status})."
            )

    return downloaded


def debug_table_structure(pdb_id):
    """
    Debug function to show the actual table structure for a PDB ID.
    """
    print(f"\n🔍 Debugging table structure for {pdb_id}")
    print("="*50)
    
    url = f"{BASE_STRUCTURE_URL}{pdb_id}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    
    try:
        resp = requests.get(url, timeout=15, headers=headers)
        resp.raise_for_status()
        
        soup = BeautifulSoup(resp.text, 'html.parser')
        tables = soup.find_all('table')
        
        print(f"Found {len(tables)} tables on the page")
        
        for table_idx, table in enumerate(tables):
            headers = table.find_all('th')
            header_text = ' '.join([th.get_text().strip() for th in headers])
            
            if any(keyword in header_text.lower() for keyword in ['molecule', 'organism', 'chains']):
                print(f"\nTable {table_idx} (potential macromolecules table):")
                print(f"Headers: {header_text}")
                
                rows = table.find_all('tr')
                print(f"Number of rows: {len(rows)}")
                
                for row_idx, row in enumerate(rows[:5]):  # Show first 5 rows
                    cells = row.find_all(['td', 'th'])
                    if cells:
                        row_text = ' | '.join([cell.get_text().strip() for cell in cells])
                        print(f"  Row {row_idx}: {row_text}")
                        
                        # Check if this row contains VDR
                        row_text_lower = row_text.lower()
                        vdr_found = any(vdr_id.lower() in row_text_lower for vdr_id in VDR_IDENTIFIERS)
                        if vdr_found:
                            print(f"    ✓ VDR FOUND in this row!")
                            print(f"    Cell details:")
                            for cell_idx, cell in enumerate(cells):
                                cell_text = cell.get_text().strip()
                                print(f"      Cell {cell_idx}: '{cell_text}'")
                
                if len(rows) > 5:
                    print(f"  ... and {len(rows)-5} more rows")
    
    except Exception as e:
        print(f"❌ Error: {e}")


def test_single_structure(pdb_id):
    """
    Test VDR organism extraction for a single PDB ID with detailed debugging.
    """
    print(f"\n🔍 Testing VDR organism extraction for {pdb_id}")
    print("="*50)
    
    # First show the table structure
    debug_table_structure(pdb_id)
    
    # Enable debug logging temporarily
    original_level = logging.getLogger().level
    logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Test the full extraction
        print(f"\n🧬 Testing VDR-specific extraction...")
        structure_info = get_complete_structure_info(pdb_id)
        
        print(f"\n📋 Final Results for {pdb_id}:")
        print(f"Structure Name: {structure_info['structure_name']}")
        print(f"VDR Organism: {structure_info['organism']}")
        print(f"Resolution: {structure_info['resolution']} Å")
        
    except Exception as e:
        print(f"❌ Error testing {pdb_id}: {e}")
    
    finally:
        # Restore original logging level
        logging.getLogger().setLevel(original_level)


def test_specific_structures():
    """
    Test VDR organism extraction for structures that should contain VDR from specific organisms.
    """
    test_structures = ["7QPI", "7OXZ", "7OY4", "8PWD", "4FHH"]
    
    print(f"\n🧪 Testing specific structures for VDR organism extraction")
    print("="*60)
    
    # Enable debug logging temporarily
    original_level = logging.getLogger().level
    logging.getLogger().setLevel(logging.DEBUG)
    
    for pdb_id in test_structures:
        print(f"\n🔍 Testing {pdb_id}:")
        print("-" * 30)
        
        try:
            name, org = scrape_structure_info(pdb_id)
            print(f"Structure name: {name}")
            print(f"VDR organism: {org}")
            
            # Also test resolution
            res = get_resolution_from_api(pdb_id)
            print(f"Resolution: {res} Å")
            
        except Exception as e:
            print(f"❌ Error: {e}")
    
    # Restore logging level
    logging.getLogger().setLevel(original_level)


def validate_csv_organisms(csv_path=CSV_FILENAME):
    """
    Validate the VDR organism distribution in the saved CSV file.
    """
    if not os.path.exists(csv_path):
        print(f"❌ CSV file {csv_path} not found")
        return
    
    print(f"\n📊 Validating VDR organism distribution in {csv_path}")
    print("="*50)
    
    organism_counts = {}
    uncertain_structures = []
    
    try:
        with open(csv_path, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            for row in reader:
                pdb_id = row['pdb_id']
                organism = row['organism']
                
                if not organism:
                    organism = 'Unknown'
                
                # Count organisms
                organism_counts[organism] = organism_counts.get(organism, 0) + 1
                
                # Track uncertain structures
                if '(uncertain)' in organism:
                    uncertain_structures.append(pdb_id)
        
        # Display results
        print(f"\n🔢 VDR organism counts in CSV:")
        for org, count in sorted(organism_counts.items(), key=lambda x: x[1], reverse=True):
            print(f"  {org}: {count} structures")
        
        if uncertain_structures:
            print(f"\n⚠️  Uncertain structures ({len(uncertain_structures)}):")
            print("  These used fallback organism detection method")
            if len(uncertain_structures) <= 10:
                print(f"  IDs: {', '.join(uncertain_structures)}")
            else:
                print(f"  First 10 IDs: {', '.join(uncertain_structures[:10])}")
        
        # Check for potential issues
        total_certain = sum(count for org, count in organism_counts.items() if '(uncertain)' not in org and org != 'Unknown')
        total_uncertain = len(uncertain_structures)
        total_unknown = organism_counts.get('Unknown', 0)
        
        print(f"\n📊 Summary:")
        print(f"  VDR-specific detection: {total_certain} structures")
        print(f"  Fallback detection: {total_uncertain} structures")
        print(f"  No organism found: {total_unknown} structures")
        
        if total_uncertain > total_certain * 0.3:  # More than 30% uncertain
            print(f"\n⚠️  WARNING: High proportion of uncertain detections ({total_uncertain}/{total_certain + total_uncertain})")
            print("   The VDR-specific detection may need improvement.")
    
    except Exception as e:
        print(f"❌ Error reading CSV: {e}")


# Updated table creation functions (simplified since no multi-organism handling)
def create_vdr_structure_table(csv_filename="vdr_structures.csv"):
    """
    Create a table showing VDR structure distribution by species and resolution.
    Uses VDR-specific organism identification.
    """
    
    # Species mapping for display names
    species_display_names = {
        "Homo sapiens": "Human",
        "Rattus norvegicus": "Rats", 
        "Danio rerio": "Zebra Fish",
        "Petromyzon marinus": "Sea Lamprey"
    }
    
    # Resolution bins
    resolution_bins = RESOLUTION_BINS
    
    # Initialize the table structure
    table_data = defaultdict(lambda: defaultdict(int))
    
    # Check if CSV file exists
    if not os.path.exists(csv_filename):
        print(f"❌ Error: {csv_filename} not found!")
        print("Please run the VDR extraction script first to generate the CSV file.")
        return
    
    # Read the CSV file
    try:
        with open(csv_filename, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            total_structures = 0
            processed_structures = 0
            
            for row in reader:
                total_structures += 1
                pdb_id = row['pdb_id']
                organism = row['organism'].strip()
                resolution_str = row['resolution'].strip()
                
                # Skip structures without resolution data
                if not resolution_str or resolution_str == '':
                    print(f"⚠️  Skipping {pdb_id}: No resolution data")
                    continue
                
                try:
                    resolution = float(resolution_str)
                except ValueError:
                    print(f"⚠️  Skipping {pdb_id}: Invalid resolution '{resolution_str}'")
                    continue
                
                # Clean organism name (remove uncertainty markers)
                clean_organism = organism.replace(" (uncertain)", "")
                
                # Determine resolution bin
                resolution_bin = None
                for bin_name, test_fn in resolution_bins.items():
                    if test_fn(resolution):
                        resolution_bin = bin_name
                        break
                
                if resolution_bin is None:
                    print(f"⚠️  {pdb_id}: Resolution {resolution} Å doesn't fit any bin")
                    continue
                
                # Count this structure if organism is in our species of interest
                if clean_organism in species_display_names:
                    display_name = species_display_names[clean_organism]
                    table_data[resolution_bin][display_name] += 1
                    processed_structures += 1
            
            print(f"📊 Processed {processed_structures} out of {total_structures} structures from {csv_filename}")
    
    except Exception as e:
        print(f"❌ Error reading {csv_filename}: {e}")
        return
    
    # Create and display the table
    print_vdr_table(table_data, species_display_names)


def print_vdr_table(table_data, species_display_names):
    """
    Print the VDR structure table in the specified format.
    """
    
    # Define the order of columns and rows
    species_order = ["Human", "Rats", "Zebra Fish", "Sea Lamprey"]
    resolution_order = ["1-2 Å", "2-2.5 Å", "2.5+ Å"]
    
    # Calculate column widths
    max_species_width = max(len(species) for species in species_order)
    col_width = max(max_species_width, 8)  # Minimum width of 8
    
    # Print header
    print(f"\n{'':>12}", end="")
    for species in species_order:
        print(f"{species:>{col_width}}", end="")
    print()  # New line
    
    # Print separator row
    print("-" * (12 + col_width * len(species_order)))
    
    # Print data rows
    totals = defaultdict(int)
    for res_bin in resolution_order:
        print(f"{res_bin:>12}", end="")
        for species in species_order:
            count = table_data[res_bin][species]
            print(f"{count if count > 0 else '-':>{col_width}}", end="")
            totals[species] += count
        print()  # New line
    
    # Print separator before totals
    print("-" * (12 + col_width * len(species_order)))
    
    # Print totals row
    print(f"{'Total:':>12}", end="")
    grand_total = 0
    for species in species_order:
        total = totals[species]
        print(f"{total if total > 0 else '-':>{col_width}}", end="")
        grand_total += total
    print()  # New line
    
    print(f"\n📈 Grand Total: {grand_total} VDR structures")
    
    # Print statistics
    print(f"\n📊 VDR Statistics:")
    for species in species_order:
        if totals[species] > 0:
            percentage = (totals[species] / grand_total) * 100
            print(f"  VDR from {species}: {totals[species]} structures ({percentage:.1f}%)")


# -----------
# Main Routine
# -----------
def main():
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

    print("🔍 VDR Structure Extraction Tool (VDR-Specific Organism Detection)")
    print("="*70)
    print("This tool identifies organisms based on the VDR protein specifically,")
    print("not other proteins that might be present in the structure.\n")
                
    # Ask what the user wants to do
    print("🎯 What would you like to do?")
    print("1. Test specific structures for VDR organism extraction")
    print("2. Test a single PDB ID")  
    print("3. Run full processing")
    
    choice = input("\nEnter choice (1-3) [default: 3]: ").strip()
    
    if choice == "1":
        test_specific_structures()
        
        continue_choice = input(f"\nContinue with full processing? (y/N): ").strip().lower()
        if continue_choice not in ['y', 'yes']:
            print("👋 Exiting. Run again when ready!")
            return
            
    elif choice == "2":
        test_pdb = input("Enter PDB ID to test (e.g., 7QPI): ").strip().upper()
        if test_pdb:
            test_single_structure(test_pdb)
            
            continue_choice = input(f"\nContinue with full processing? (y/N): ").strip().lower()
            if continue_choice not in ['y', 'yes']:
                print("👋 Exiting. Run again when ready!")
                return

    print(f"\n📋 Species detection priority: {' > '.join(SPECIES_LIST)}")
    print("🎯 Focus: Detecting organisms of VDR protein specifically")

    # Step 1: Query the PDB
    logging.info("Querying PDB database...")
    all_results = query_pdb_structures(query)
    if not all_results:
        logging.info("No results returned from PDB query. Exiting.")
        return

    # Step 2: Filter out DNA-binding entries
    filtered = filter_pdb_results(all_results, NOT_ALLOWED_PDB_IDS)
    if not filtered:
        logging.info("No structures found after filtering. Exiting.")
        return

    # Step 3: Print filtered results and collect IDs
    pdb_ids_to_process = []
    print(f"\n📋 Found {len(filtered)} structures (filtered):")
    for entry in filtered:
        pdb_id = entry["pdb_id"]
        score = entry["score"]
        print(f"  • PDB ID: {pdb_id}, Score: {score}")
        pdb_ids_to_process.append(pdb_id)

    print(f"\n📊 Processing {len(pdb_ids_to_process)} structures...")
    print(f"Rate limiting: {REQUEST_DELAY} seconds between requests")
    print(f"Maximum retries per request: {MAX_RETRIES}")

    # Step 4: Process structures and save to CSV
    all_structures = process_and_save_structures(pdb_ids_to_process)

    # Step 4.5: Validate the CSV results
    validate_csv_organisms()

    # Step 5: Download structure files (optional)
    download_choice = input(f"\n❓ Download structure files to '{OUTPUT_DIR}'? (y/N): ").strip().lower()
    if download_choice in ['y', 'yes']:
        print(f"\n📥 Downloading structure files...")
        downloaded_info = download_structures(pdb_ids_to_process, OUTPUT_DIR)

        if downloaded_info:
            print(f"\n✅ Successfully downloaded {len(downloaded_info)} files:")
            for pid, ext in downloaded_info:
                print(f"  ✔ {pid}{ext}")
        else:
            print("\n❌ No files were successfully downloaded.")

        # Check for any IDs that didn't succeed
        downloaded_ids = {pid for pid, _ in downloaded_info}
        missing = set(pdb_ids_to_process) - downloaded_ids
        if missing:
            print(f"\n⚠️  Failed to download {len(missing)} files:")
            for mid in missing:
                print(f"  ✘ {mid}")

    create_table = input(f"\n❓ Create VDR structure distribution table? (y/N): ").strip().lower()
    if create_table not in ['n', 'no']:
        print(f"\n📊 Creating VDR structure distribution table...")
        create_vdr_structure_table(CSV_FILENAME)
            
    print(f"\n🎉 Processing complete! Check '{CSV_FILENAME}' for VDR structure data.")
    print("📋 Organisms now represent the species of the VDR protein specifically.")

if __name__ == "__main__":
    main()