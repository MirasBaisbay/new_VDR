
# PyMOL script for analyzing repaired 5H1E structure
load alphafold_repair/5H1E_alphafold_repaired.pdb, repaired_5h1e

# Color by source
color lightblue, repaired_5h1e and resi 123-159    # Experimental
color orange, repaired_5h1e and resi 160-217       # AlphaFold grafted
color lightblue, repaired_5h1e and resi 218-419    # Experimental

# Highlight problem areas
color red, repaired_5h1e and resi 159-160          # Junction 1
color red, repaired_5h1e and resi 217-218          # Junction 2

# Show cartoon representation
hide everything
show cartoon, repaired_5h1e

# Add distance measurements
distance dist1, repaired_5h1e and resi 159 and name CA, repaired_5h1e and resi 160 and name CA
distance dist2, repaired_5h1e and resi 217 and name CA, repaired_5h1e and resi 218 and name CA

# Center on junction regions
center repaired_5h1e and resi 155-165
zoom repaired_5h1e and resi 155-165

# Save image
png junction1_analysis.png, dpi=300

center repaired_5h1e and resi 213-223
zoom repaired_5h1e and resi 213-223
png junction2_analysis.png, dpi=300

# Full structure view
zoom repaired_5h1e
png full_structure.png, dpi=300
