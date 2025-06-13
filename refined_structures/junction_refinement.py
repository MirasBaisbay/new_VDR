
# MODELLER script for junction refinement
from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = ['.']

# Refine only junction regions
class JunctionRefine(loopmodel):
    def select_loop_atoms(self):
        # Select junction regions for refinement
        return (selection(self.residue_range('155:A', '165:A')) |
                selection(self.residue_range('213:A', '223:A')))

m = JunctionRefine(env,
                  alnfile='junction_refine.ali',
                  knowns=('template'),
                  sequence='5h1e_repaired')

m.loop.starting_model = 1
m.loop.ending_model = 5
m.loop.md_level = refine.fast

m.make()
