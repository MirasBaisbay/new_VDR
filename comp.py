#!/usr/bin/env python3
"""
Batch Structure Analyzer - Check all 48 human VDR structures
"""

import os
import glob
from pdbfixer import PDBFixer

class BatchStructureAnalyzer:
    def __init__(self, original_dir="pdb_files/Homo sapiens", repaired_dir="vdr_processed/repaired_structures"):
        self.original_dir = original_dir
        self.repaired_dir = repaired_dir
        
        self.aa_codes = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
    
    def analyze_single_structure(self, pdb_file):
        """
        Analyze a single structure quickly
        """
        try:
            fixer = PDBFixer(filename=pdb_file)
            
            total_atoms = 0
            total_residues = 0
            protein_residues = 0
            hetero_residues = 0
            water_residues = 0
            chain_count = 0
            
            heterogens = set()
            protein_sequence = ""
            
            for chain in fixer.topology.chains():
                chain_count += 1
                
                for residue in chain.residues():
                    total_residues += 1
                    total_atoms += len(list(residue.atoms()))
                    
                    if residue.name in self.aa_codes:
                        protein_residues += 1
                        protein_sequence += self.aa_codes[residue.name]
                    elif residue.name in ['HOH', 'WAT']:
                        water_residues += 1
                    else:
                        hetero_residues += 1
                        heterogens.add(residue.name)
            
            return {
                'success': True,
                'total_atoms': total_atoms,
                'total_residues': total_residues,
                'protein_residues': protein_residues,
                'hetero_residues': hetero_residues,
                'water_residues': water_residues,
                'chain_count': chain_count,
                'protein_sequence': protein_sequence,
                'heterogens': sorted(list(heterogens))
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }
    
    def batch_analyze_all_structures(self):
        """
        Analyze all structures in batch
        """
        print("BATCH STRUCTURE ANALYSIS - ALL 48 HUMAN VDR STRUCTURES")
        print("="*70)
        
        # Find all original structures
        original_files = glob.glob(os.path.join(self.original_dir, "*.pdb"))
        pdb_codes = [os.path.basename(f).replace('.pdb', '') for f in original_files]
        
        if not pdb_codes:
            print(f"No PDB files found in {self.original_dir}")
            return
        
        print(f"Found {len(pdb_codes)} structures to analyze")
        print("")
        
        results = {}
        successful_comparisons = 0
        failed_comparisons = 0
        
        # Analyze each structure
        for i, pdb_code in enumerate(sorted(pdb_codes), 1):
            print(f"[{i:2d}/{len(pdb_codes)}] Analyzing {pdb_code}...", end=" ")
            
            original_file = os.path.join(self.original_dir, f"{pdb_code}.pdb")
            repaired_file = os.path.join(self.repaired_dir, f"{pdb_code}_repaired.pdb")
            
            if not os.path.exists(repaired_file):
                print("❌ Repaired file not found")
                results[pdb_code] = {'status': 'no_repaired_file'}
                failed_comparisons += 1
                continue
            
            # Analyze both structures
            orig_result = self.analyze_single_structure(original_file)
            rep_result = self.analyze_single_structure(repaired_file)
            
            if not orig_result['success']:
                print(f"❌ Original analysis failed: {orig_result['error']}")
                failed_comparisons += 1
                continue
            
            if not rep_result['success']:
                print(f"❌ Repaired analysis failed: {rep_result['error']}")
                failed_comparisons += 1
                continue
            
            # Compare results
            comparison = self.compare_structures(orig_result, rep_result)
            results[pdb_code] = {
                'status': 'success',
                'original': orig_result,
                'repaired': rep_result,
                'comparison': comparison
            }
            
            # Quick status
            if comparison['protein_sequence_identical']:
                if comparison['atoms_added'] > 0:
                    print(f"✅ Perfect (+{comparison['atoms_added']} atoms)")
                else:
                    print(f"✅ Identical")
            else:
                print(f"⚠️  Sequence changed ({comparison['protein_residue_diff']:+d} residues)")
            
            successful_comparisons += 1
        
        # Generate summary report
        self.generate_batch_summary(results, successful_comparisons, failed_comparisons)
    
    def compare_structures(self, orig_result, rep_result):
        """
        Compare original vs repaired structure
        """
        return {
            'atoms_added': rep_result['total_atoms'] - orig_result['total_atoms'],
            'total_residue_diff': rep_result['total_residues'] - orig_result['total_residues'],
            'protein_residue_diff': rep_result['protein_residues'] - orig_result['protein_residues'],
            'protein_sequence_identical': orig_result['protein_sequence'] == rep_result['protein_sequence'],
            'chain_count_diff': rep_result['chain_count'] - orig_result['chain_count'],
            'heterogens_original': orig_result['heterogens'],
            'heterogens_repaired': rep_result['heterogens']
        }
    
    def generate_batch_summary(self, results, successful, failed):
        """
        Generate comprehensive summary report
        """
        print(f"\n" + "="*70)
        print("BATCH ANALYSIS SUMMARY")
        print("="*70)
        
        total_structures = len(results)
        print(f"Total structures: {total_structures}")
        print(f"Successful comparisons: {successful}")
        print(f"Failed comparisons: {failed}")
        
        if successful == 0:
            print("No successful comparisons to analyze!")
            return
        
        # Statistics
        successful_results = [r for r in results.values() if r['status'] == 'success']
        
        # Atom addition statistics
        atoms_added = [r['comparison']['atoms_added'] for r in successful_results]
        min_atoms = min(atoms_added)
        max_atoms = max(atoms_added)
        avg_atoms = sum(atoms_added) / len(atoms_added)
        
        print(f"\nATOM ADDITION STATISTICS:")
        print(f"  Minimum atoms added: {min_atoms}")
        print(f"  Maximum atoms added: {max_atoms}")
        print(f"  Average atoms added: {avg_atoms:.0f}")
        
        # Sequence preservation check
        identical_sequences = sum(1 for r in successful_results if r['comparison']['protein_sequence_identical'])
        print(f"\nSEQUENCE PRESERVATION:")
        print(f"  Identical protein sequences: {identical_sequences}/{successful}")
        print(f"  Sequence preservation rate: {100*identical_sequences/successful:.1f}%")
        
        if identical_sequences < successful:
            changed_sequences = [code for code, r in results.items() 
                               if r['status'] == 'success' and not r['comparison']['protein_sequence_identical']]
            print(f"  Structures with changed sequences: {changed_sequences}")
        
        # Chain organization
        chain_changes = [r['comparison']['chain_count_diff'] for r in successful_results]
        reorganized = sum(1 for diff in chain_changes if diff != 0)
        print(f"\nCHAIN ORGANIZATION:")
        print(f"  Structures with chain reorganization: {reorganized}")
        
        # Heterogen analysis
        all_heterogens = set()
        for r in successful_results:
            all_heterogens.update(r['comparison']['heterogens_original'])
            all_heterogens.update(r['comparison']['heterogens_repaired'])
        
        print(f"\nHETEROGEN ANALYSIS:")
        print(f"  Unique heterogens found: {sorted(all_heterogens)}")
        
        # Detailed results table
        print(f"\nDETAILED RESULTS:")
        print(f"{'PDB':<6} {'Original':<8} {'Repaired':<8} {'+Atoms':<7} {'Chains':<7} {'SeqOK':<6} {'Heterogens'}")
        print("-" * 70)
        
        for pdb_code in sorted(results.keys()):
            result = results[pdb_code]
            
            if result['status'] != 'success':
                print(f"{pdb_code:<6} {'FAILED':<8} {'FAILED':<8} {'N/A':<7} {'N/A':<7} {'N/A':<6} N/A")
                continue
            
            orig = result['original']
            rep = result['repaired']
            comp = result['comparison']
            
            seq_ok = "✓" if comp['protein_sequence_identical'] else "✗"
            chains = f"{orig['chain_count']}→{rep['chain_count']}"
            heterogens = ','.join(comp['heterogens_original'][:3])  # First 3 only
            
            print(f"{pdb_code:<6} {orig['protein_residues']:<8} {rep['protein_residues']:<8} "
                  f"{comp['atoms_added']:<7} {chains:<7} {seq_ok:<6} {heterogens}")
        
        # Save detailed report
        self.save_detailed_report(results)
    
    def save_detailed_report(self, results):
        """
        Save detailed report to file
        """
        report_file = "batch_structure_analysis_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("Batch Structure Analysis Report\n")
            f.write("="*50 + "\n\n")
            
            f.write("SUMMARY:\n")
            f.write(f"Total structures analyzed: {len(results)}\n")
            
            successful_results = [r for r in results.values() if r['status'] == 'success']
            f.write(f"Successful comparisons: {len(successful_results)}\n\n")
            
            f.write("DETAILED RESULTS:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'PDB Code':<10} {'Orig Atoms':<12} {'Rep Atoms':<12} {'Atoms Added':<12} "
                   f"{'Prot Res':<10} {'Seq OK':<8} {'Heterogens'}\n")
            f.write("-" * 80 + "\n")
            
            for pdb_code in sorted(results.keys()):
                result = results[pdb_code]
                
                if result['status'] != 'success':
                    f.write(f"{pdb_code:<10} {'FAILED':<12} {'FAILED':<12} {'N/A':<12} "
                           f"{'N/A':<10} {'N/A':<8} N/A\n")
                    continue
                
                orig = result['original']
                rep = result['repaired']
                comp = result['comparison']
                
                seq_ok = "YES" if comp['protein_sequence_identical'] else "NO"
                heterogens = ','.join(comp['heterogens_original'])
                
                f.write(f"{pdb_code:<10} {orig['total_atoms']:<12} {rep['total_atoms']:<12} "
                       f"{comp['atoms_added']:<12} {orig['protein_residues']:<10} "
                       f"{seq_ok:<8} {heterogens}\n")
            
            # Problem structures
            problem_structures = [code for code, r in results.items() 
                                if r['status'] == 'success' and not r['comparison']['protein_sequence_identical']]
            
            if problem_structures:
                f.write(f"\nSTRUCTURES WITH SEQUENCE CHANGES:\n")
                for code in problem_structures:
                    result = results[code]
                    comp = result['comparison']
                    f.write(f"{code}: {comp['protein_residue_diff']:+d} protein residues\n")
        
        print(f"\nDetailed report saved: {report_file}")

def main():
    analyzer = BatchStructureAnalyzer()
    analyzer.batch_analyze_all_structures()

if __name__ == "__main__":
    main()