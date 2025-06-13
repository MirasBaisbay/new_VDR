#!/usr/bin/env python3
"""
Simple Ramachandran Plotter Runner

Just loops through PDB files and runs RamachandranPlotter.py on each one.
"""

import os
import subprocess
import glob
from pathlib import Path

# Configuration - EDIT THESE PATHS
pdb_directory = "pdb_files"  # Directory containing your PDB files
output_directory = "ramachandran_results"  # Where to save results
ramachandran_script = r"C:\Users\Meiras\Desktop\Research\Research Ferdinand\new_VDR\Ramachandran_Plotter-main\RamachandranPlotter.py"

# Create output directory
os.makedirs(output_directory, exist_ok=True)

# Find all PDB files
pdb_files = []
for pattern in ["*.pdb", "*.PDB"]:
    pdb_files.extend(glob.glob(os.path.join(pdb_directory, pattern)))

print(f"Found {len(pdb_files)} PDB files")

# Process each PDB file
successful = 0
failed = 0

for i, pdb_file in enumerate(pdb_files, 1):
    pdb_name = os.path.basename(pdb_file)
    print(f"\n[{i}/{len(pdb_files)}] Processing: {pdb_name}")
    
    try:
        # Run RamachandranPlotter.py
        cmd = [
            "python", 
            ramachandran_script,
            "--pdb", pdb_file,
            "--out_dir", output_directory,
            "--plot_type", "0",  # All residues
            "--save_csv"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print(f"  ✓ Success")
            successful += 1
        else:
            print(f"  ✗ Failed - Return code: {result.returncode}")
            if result.stderr:
                print(f"    Error: {result.stderr[:100]}")
            failed += 1
            
    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout (>60 seconds)")
        failed += 1
    except Exception as e:
        print(f"  ✗ Error: {e}")
        failed += 1

print(f"\n" + "="*50)
print(f"SUMMARY:")
print(f"Successful: {successful}")
print(f"Failed: {failed}")
print(f"Total: {len(pdb_files)}")
print(f"Results saved to: {output_directory}")

# Check what files were created
plot_files = glob.glob(os.path.join(output_directory, "*.png"))
csv_files = glob.glob(os.path.join(output_directory, "*.csv"))

print(f"\nOutput files created:")
print(f"Plots: {len(plot_files)}")
print(f"CSV files: {len(csv_files)}")