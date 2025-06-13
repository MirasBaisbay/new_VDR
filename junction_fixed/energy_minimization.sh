
# GROMACS energy minimization script for junction repair
# Run after junction correction to relax local strain

# 1. Convert PDB to GROMACS format
gmx pdb2gmx -f 5H1E_junction_fixed_v2.pdb -o 5h1e_processed.gro -p 5h1e.top -ff amber99sb-ildn -water tip3p

# 2. Create simulation box
gmx editconf -f 5h1e_processed.gro -o 5h1e_box.gro -c -d 1.0 -bt cubic

# 3. Add solvent
gmx solvate -cp 5h1e_box.gro -cs spc216.gro -o 5h1e_solv.gro -p 5h1e.top

# 4. Energy minimization (steepest descent)
gmx grompp -f em.mdp -c 5h1e_solv.gro -p 5h1e.top -o em.tpr
gmx mdrun -v -deffnm em

# 5. Extract minimized structure
gmx editconf -f em.gro -o 5h1e_minimized.pdb
