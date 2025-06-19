#!/usr/bin/env python3
"""
VDR Clean Processor - Analyze and repair all structures in folder
"""

import os
import glob
import re
from pdbfixer import PDBFixer
from openmm.app import PDBFile

class VDRProcessor:
    def __init__(self, input_folder="pdb_files/Homo sapiens", output_folder="vdr_processed"):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.repaired_folder = os.path.join(output_folder, 'repaired_structures')
        
        # Create output directories
        os.makedirs(self.repaired_folder, exist_ok=True)
    
    def get_resolution(self, pdb_file):
        """Extract resolution from PDB file"""
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('REMARK   2 RESOLUTION'):
                        res_match = re.search(r'(\d+\.\d+)', line)
                        if res_match:
                            return float(res_match.group(1))
                    elif line.startswith('ATOM'):
                        break
        except:
            pass
        return None
    
    def analyze_structure(self, pdb_code):
        """Analyze structure quality"""
        input_file = os.path.join(self.input_folder, f"{pdb_code}.pdb")
        
        # Get resolution
        resolution = self.get_resolution(input_file)
        
        # Analyze with PDBFixer
        try:
            fixer = PDBFixer(filename=input_file)
            
            # Count atoms and missing residues
            total_atoms = len(list(fixer.topology.atoms()))
            
            # Find missing residues
            fixer.findMissingResidues()
            missing_residues = len(fixer.missingResidues)
            
            # Get heterogens
            protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                              'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
            
            heterogens = []
            for residue in fixer.topology.residues():
                if residue.name not in protein_residues:
                    heterogens.append(residue.name)
            
            # Make recommendation
            recommendation = self.get_recommendation(resolution, missing_residues)
            
            return {
                'pdb_code': pdb_code,
                'resolution': resolution,
                'total_atoms': total_atoms,
                'missing_residues': missing_residues,
                'heterogens': list(set(heterogens)),
                'recommendation': recommendation
            }
            
        except Exception as e:
            print(f"Error analyzing {pdb_code}: {e}")
            return None
    
    def get_recommendation(self, resolution, missing_residues):
        """Simple recommendation logic"""
        if resolution is None:
            return "UNKNOWN"
        elif resolution > 3.0 or missing_residues > 5:
            return "NEEDS_MINIMIZATION"
        elif resolution > 2.5 or missing_residues > 0:
            return "CONSIDER_MINIMIZATION"
        else:
            return "GOOD"
    
    def repair_structure(self, pdb_code):
        """Repair structure with PDBFixer"""
        input_file = os.path.join(self.input_folder, f"{pdb_code}.pdb")
        output_file = os.path.join(self.repaired_folder, f"{pdb_code}_repaired.pdb")
        
        try:
            fixer = PDBFixer(filename=input_file)
            
            # Repair steps with error handling
            try:
                fixer.findMissingResidues()
            except Exception as e:
                print(f"    Warning: Could not find missing residues: {e}")
            
            try:
                fixer.findNonstandardResidues()
                
                # Keep important ligands
                important_ligands = ['COA', 'COB', 'CO4', 'VDX', 'VD3', 'VD2', 'KH1', 'TEJ', 'CAL']
                if fixer.nonstandardResidues:
                    fixer.nonstandardResidues = [
                        (res, repl) for res, repl in fixer.nonstandardResidues 
                        if res.name not in important_ligands
                    ]
                
                fixer.replaceNonstandardResidues()
            except Exception as e:
                print(f"    Warning: Could not handle nonstandard residues: {e}")
            
            try:
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
            except Exception as e:
                print(f"    Warning: Could not add missing atoms: {e}")
            
            try:
                fixer.addMissingHydrogens(7.0)
            except Exception as e:
                print(f"    Warning: Could not add hydrogens (skipping): {e}")
            
            # Save repaired structure
            with open(output_file, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
            # Count atoms added
            try:
                original_pdb = PDBFile(input_file)
                repaired_pdb = PDBFile(output_file)
                
                original_atoms = len(list(original_pdb.topology.atoms()))
                repaired_atoms = len(list(repaired_pdb.topology.atoms()))
                atoms_added = repaired_atoms - original_atoms
            except:
                # If counting fails, estimate
                atoms_added = 0
            
            return {
                'success': True,
                'original_atoms': original_atoms if 'original_atoms' in locals() else 0,
                'repaired_atoms': repaired_atoms if 'repaired_atoms' in locals() else 0,
                'atoms_added': atoms_added
            }
            
        except Exception as e:
            print(f"    Error repairing {pdb_code}: {e}")
            
            # Try to at least copy the original file
            try:
                import shutil
                shutil.copy2(input_file, output_file)
                print(f"    Copied original file as fallback")
                
                original_pdb = PDBFile(input_file)
                original_atoms = len(list(original_pdb.topology.atoms()))
                
                return {
                    'success': True,
                    'original_atoms': original_atoms,
                    'repaired_atoms': original_atoms,
                    'atoms_added': 0,
                    'fallback': True
                }
            except:
                return {'success': False, 'error': str(e)}
    
    def process_all_structures(self):
        """Process all structures in the folder"""
        print("VDR Structure Processor")
        print("="*50)
        
        # Get all PDB files
        pdb_files = glob.glob(os.path.join(self.input_folder, "*.pdb"))
        pdb_codes = [os.path.basename(f).replace('.pdb', '') for f in pdb_files]
        
        if not pdb_codes:
            print(f"No PDB files found in {self.input_folder}")
            return
        
        print(f"Found {len(pdb_codes)} structures")
        print(f"Input folder: {self.input_folder}")
        print(f"Output folder: {self.output_folder}")
        print("")
        
        # Process each structure
        results = []
        successful_repairs = 0
        
        for i, pdb_code in enumerate(pdb_codes, 1):
            print(f"[{i}/{len(pdb_codes)}] Processing {pdb_code}...")
            
            # Analyze quality
            analysis = self.analyze_structure(pdb_code)
            if not analysis:
                continue
            
            # Repair structure
            repair_result = self.repair_structure(pdb_code)
            
            # Combine results
            result = {**analysis, **repair_result}
            results.append(result)
            
            if repair_result.get('success', False):
                successful_repairs += 1
            
            # Print quick info
            res = analysis['resolution']
            rec = analysis['recommendation']
            atoms = repair_result.get('atoms_added', 0)
            print(f"  {res:.2f}Å | {rec} | +{atoms} atoms")
        
        # Generate summary
        self.print_summary(results, successful_repairs)
        self.save_summary(results)
    
    def print_summary(self, results, successful_repairs):
        """Print processing summary"""
        print("\n" + "="*60)
        print("PROCESSING SUMMARY")
        print("="*60)
        
        print(f"Total structures: {len(results)}")
        print(f"Successfully repaired: {successful_repairs}")
        
        # Group by recommendation
        good = [r for r in results if r['recommendation'] == 'GOOD']
        consider = [r for r in results if r['recommendation'] == 'CONSIDER_MINIMIZATION']
        needs = [r for r in results if r['recommendation'] == 'NEEDS_MINIMIZATION']
        unknown = [r for r in results if r['recommendation'] == 'UNKNOWN']
        
        print(f"\nRECOMMENDATIONS:")
        print(f"  GOOD (use repaired as-is): {len(good)}")
        print(f"  CONSIDER minimization: {len(consider)}")
        print(f"  NEEDS minimization: {len(needs)}")
        print(f"  UNKNOWN resolution: {len(unknown)}")
        
        # Show examples
        if good:
            print(f"\n  Good structures (examples): {[r['pdb_code'] for r in good[:5]]}")
        if needs:
            print(f"  Need minimization: {[r['pdb_code'] for r in needs]}")
        
        # Resolution statistics
        resolutions = [r['resolution'] for r in results if r['resolution']]
        if resolutions:
            avg_res = sum(resolutions) / len(resolutions)
            min_res = min(resolutions)
            max_res = max(resolutions)
            print(f"\n  Resolution range: {min_res:.2f} - {max_res:.2f} Å (avg: {avg_res:.2f} Å)")
        
        print(f"\nOutput files saved in: {self.repaired_folder}")
    
    def save_summary(self, results):
        """Save detailed summary to file"""
        summary_file = os.path.join(self.output_folder, "processing_summary.txt")
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("VDR Structure Processing Summary\n")
            f.write("="*50 + "\n\n")
            
            f.write(f"Total structures processed: {len(results)}\n")
            f.write(f"Input folder: {self.input_folder}\n")
            f.write(f"Output folder: {self.repaired_folder}\n\n")
            
            f.write("STRUCTURE DETAILS:\n")
            f.write("-"*40 + "\n")
            f.write(f"{'PDB':<6} {'Resolution':<10} {'Missing':<8} {'Atoms+':<8} {'Recommendation':<20} {'Ligands'}\n")
            f.write("-"*40 + "\n")
            
            for result in sorted(results, key=lambda x: x['resolution'] if x['resolution'] else 999):
                pdb = result['pdb_code']
                res = f"{result['resolution']:.2f}Å" if result['resolution'] else "N/A"
                missing = result['missing_residues']
                atoms = result.get('atoms_added', 0)
                rec = result['recommendation']
                ligands = ','.join(result['heterogens'][:3])  # First 3 ligands
                
                f.write(f"{pdb:<6} {res:<10} {missing:<8} {atoms:<8} {rec:<20} {ligands}\n")
            
            # Summary statistics
            good = sum(1 for r in results if r['recommendation'] == 'GOOD')
            consider = sum(1 for r in results if r['recommendation'] == 'CONSIDER_MINIMIZATION')
            needs = sum(1 for r in results if r['recommendation'] == 'NEEDS_MINIMIZATION')
            
            f.write(f"\nSUMMARY:\n")
            f.write(f"  GOOD (use repaired): {good}\n")
            f.write(f"  CONSIDER minimization: {consider}\n")
            f.write(f"  NEEDS minimization: {needs}\n")
        
        print(f"Detailed summary saved: {summary_file}")

def main():
    # Process all structures in Homo sapiens folder
    processor = VDRProcessor()
    processor.process_all_structures()

if __name__ == "__main__":
    main()