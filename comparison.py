import os
import statistics

def extract_bfactors_original_method(pdb_file_path):
    """Original method - only ATOM records"""
    try:
        b_factors = []
        
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):  # Only ATOM, no HETATM
                    try:
                        bfactor_str = line[60:66].strip()
                        if bfactor_str:
                            b_factor = float(bfactor_str)
                            b_factors.append(b_factor)
                    except (ValueError, IndexError):
                        continue
        
        if not b_factors:
            return None
        
        return {
            'mean_bfactor': round(statistics.mean(b_factors), 2),
            'min_bfactor': round(min(b_factors), 2),
            'max_bfactor': round(max(b_factors), 2),
            'num_atoms': len(b_factors)
        }
        
    except Exception as e:
        return None

def extract_bfactors_new_method(pdb_file_path):
    """New method - ATOM + HETATM, first model only, primary conformations"""
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
                        # Check alternative location indicator
                        altloc = line[16] if len(line) > 16 else ' '
                        if altloc in [' ', 'A']:  # Only primary or 'A' conformations
                            bfactor_str = line[60:66].strip()
                            if bfactor_str:
                                b_factor = float(bfactor_str)
                                b_factors.append(b_factor)
                    except (ValueError, IndexError):
                        continue
        
        if not b_factors:
            return None
        
        return {
            'mean_bfactor': round(statistics.mean(b_factors), 2),
            'min_bfactor': round(min(b_factors), 2),
            'max_bfactor': round(max(b_factors), 2),
            'num_atoms': len(b_factors)
        }
        
    except Exception as e:
        return None

def extract_bfactors_protein_only(pdb_file_path):
    """Protein-only method - ATOM records, first model, primary conformations"""
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
                elif in_model_1 and line.startswith('ATOM'):  # Only protein atoms
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
        
        return {
            'mean_bfactor': round(statistics.mean(b_factors), 2),
            'min_bfactor': round(min(b_factors), 2),
            'max_bfactor': round(max(b_factors), 2),
            'num_atoms': len(b_factors)
        }
        
    except Exception as e:
        return None

def analyze_pdb_content(pdb_file_path):
    """Analyze what's in a PDB file to understand differences"""
    try:
        stats = {
            'atom_count': 0,
            'hetatm_count': 0,
            'models': set(),
            'alt_conformations': set(),
            'chains': set(),
            'residue_types': set()
        }
        
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    model_num = int(line.split()[1]) if len(line.split()) > 1 else 1
                    stats['models'].add(model_num)
                elif line.startswith('ATOM'):
                    stats['atom_count'] += 1
                    if len(line) > 21:
                        stats['chains'].add(line[21])
                    if len(line) > 16:
                        altloc = line[16]
                        if altloc != ' ':
                            stats['alt_conformations'].add(altloc)
                    if len(line) > 19:
                        res_name = line[17:20].strip()
                        stats['residue_types'].add(res_name)
                elif line.startswith('HETATM'):
                    stats['hetatm_count'] += 1
                    if len(line) > 19:
                        res_name = line[17:20].strip()
                        stats['residue_types'].add(res_name)
        
        return stats
        
    except Exception as e:
        return None

def compare_bfactor_methods(pdb_file_path):
    """Compare all three B-factor extraction methods"""
    print(f"\n🔍 Analyzing B-factor extraction methods for: {os.path.basename(pdb_file_path)}")
    print("="*60)
    
    # Analyze PDB content first
    content = analyze_pdb_content(pdb_file_path)
    if content:
        print(f"📊 PDB Content Analysis:")
        print(f"  ATOM records: {content['atom_count']}")
        print(f"  HETATM records: {content['hetatm_count']}")
        print(f"  Models: {sorted(content['models']) if content['models'] else ['No MODEL records (single model)']}")
        print(f"  Alternative conformations: {sorted(content['alt_conformations']) if content['alt_conformations'] else ['None']}")
        print(f"  Chains: {sorted(content['chains'])}")
        print(f"  Unique residues: {len(content['residue_types'])} types")
    
    print(f"\n📈 B-factor Comparison:")
    print(f"{'Method':<20} | {'Mean':<8} | {'Min':<8} | {'Max':<8} | {'Atoms':<8}")
    print("-" * 60)
    
    # Method 1: Original (ATOM only)
    original = extract_bfactors_original_method(pdb_file_path)
    if original:
        print(f"{'Original (ATOM)':<20} | {original['mean_bfactor']:<8} | {original['min_bfactor']:<8} | {original['max_bfactor']:<8} | {original['num_atoms']:<8}")
    else:
        print(f"{'Original (ATOM)':<20} | {'Failed':<8} | {'Failed':<8} | {'Failed':<8} | {'Failed':<8}")
    
    # Method 2: New (ATOM + HETATM)
    new_method = extract_bfactors_new_method(pdb_file_path)
    if new_method:
        print(f"{'New (ATOM+HETATM)':<20} | {new_method['mean_bfactor']:<8} | {new_method['min_bfactor']:<8} | {new_method['max_bfactor']:<8} | {new_method['num_atoms']:<8}")
    else:
        print(f"{'New (ATOM+HETATM)':<20} | {'Failed':<8} | {'Failed':<8} | {'Failed':<8} | {'Failed':<8}")
    
    # Method 3: Protein-only (ATOM, filtered)
    protein_only = extract_bfactors_protein_only(pdb_file_path)
    if protein_only:
        print(f"{'Protein-only':<20} | {protein_only['mean_bfactor']:<8} | {protein_only['min_bfactor']:<8} | {protein_only['max_bfactor']:<8} | {protein_only['num_atoms']:<8}")
    else:
        print(f"{'Protein-only':<20} | {'Failed':<8} | {'Failed':<8} | {'Failed':<8} | {'Failed':<8}")
    
    # Analysis
    if original and new_method:
        diff = abs(original['mean_bfactor'] - new_method['mean_bfactor'])
        atom_diff = new_method['num_atoms'] - original['num_atoms']
        print(f"\n💡 Analysis:")
        print(f"  Mean B-factor difference: {diff:.2f}")
        print(f"  Additional atoms in new method: {atom_diff}")
        
        if atom_diff > 0:
            print(f"  → New method includes {atom_diff} more atoms (likely HETATM records)")
            if content and content['hetatm_count'] > 0:
                print(f"  → PDB contains {content['hetatm_count']} HETATM records (waters, ligands, etc.)")
        
        if content and content['alt_conformations']:
            print(f"  → PDB has alternative conformations: {sorted(content['alt_conformations'])}")
            print(f"  → New method filters to primary/A conformations only")
        
        if content and len(content['models']) > 1:
            print(f"  → PDB has multiple models: {sorted(content['models'])}")
            print(f"  → New method uses only first model")

def test_multiple_pdbs(pdb_directory="pdb_files"):
    """Test multiple PDB files to see patterns"""
    if not os.path.exists(pdb_directory):
        print(f"❌ Directory {pdb_directory} not found!")
        return
    
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]
    
    if not pdb_files:
        print(f"❌ No PDB files found in {pdb_directory}")
        return
    
    print(f"🔍 Testing {len(pdb_files)} PDB files...")
    
    significant_diffs = []
    
    for pdb_file in pdb_files[:5]:  # Test first 5 files
        pdb_path = os.path.join(pdb_directory, pdb_file)
        
        original = extract_bfactors_original_method(pdb_path)
        new_method = extract_bfactors_new_method(pdb_path)
        
        if original and new_method:
            diff = abs(original['mean_bfactor'] - new_method['mean_bfactor'])
            if diff > 1.0:  # Significant difference
                significant_diffs.append((pdb_file, diff, original['mean_bfactor'], new_method['mean_bfactor']))
    
    if significant_diffs:
        print(f"\n⚠️  Files with significant B-factor differences (>1.0):")
        for pdb_file, diff, orig, new in significant_diffs:
            print(f"  {pdb_file}: Original={orig}, New={new}, Diff={diff:.2f}")
        
        # Analyze the first one in detail
        first_file = os.path.join(pdb_directory, significant_diffs[0][0])
        compare_bfactor_methods(first_file)

# Usage examples:
if __name__ == "__main__":
    # Test a specific PDB file
    # compare_bfactor_methods("pdb_files/1ABC.pdb")
    
    # Test multiple files to find patterns
    test_multiple_pdbs()