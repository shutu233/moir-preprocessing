import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.io.ase import AseAtomsAdaptor

from ase.io import read, write
from ase import Atoms
import os

# import function calculating diff using StructureMatcher
import sys
sys.path.append('/global/homes/l/ltshu/Simulation/struc_compare/structurematcher')
import struc_Matcher

# Rattle method
import ase_rattle

import atoms_minimize

# unique local structure and corresponding potential energy
unique_local_struc = []
unique_local_eng = []
local_eng = []

# read atoms
at_path = "/global/homes/l/ltshu/Simulation/Structures/01-4.00_deg-540_atoms.xyz"
at = read(at_path)

N = 200
rattled_atoms = ase_rattle.Rattle(at_path, N)

# relaxation & check uniqueness of local minima
for i in range(0, N):
    rattled_atom = rattled_atoms[i]
    t = atoms_minimize.Atoms_minimize(rattled_atom)
    relaxed_rattled_atom = t[0]
    relaxed_rattled_eng = t[1]
    local_eng.append(relaxed_rattled_eng)

    if len(unique_local_struc) == 0:
        unique_local_struc.append(relaxed_rattled_atom)
        unique_local_eng.append(relaxed_rattled_eng)
    else:
        for j in range(0, len(unique_local_struc)):
            check_struc = unique_local_struc[j]
            struc_diff = struc_Matcher.structure_matcher(check_struc, relaxed_rattled_atom)
            if struc_diff[0] < 1e-2 or struc_diff[1] < 1e-2: # or 1e-3? how to choose the threshold properly
                break
            else:
                if j == len(unique_local_struc) - 1:
                    unique_local_struc.append(relaxed_rattled_atom)
                    unique_local_eng.append(relaxed_rattled_eng)
                else:
                    continue


# write the results to a file
with open("ase_rattle_MoS2WSe2_2D540atoms_4deg_N200.txt", 'w') as f:
    f.write(f"{len(unique_local_struc)}\n")
    f.write(f"{unique_local_eng}\n")


for i in range(0, len(unique_local_struc)):
    filename = f'local_{i}.xyz'
    out_dir = '/global/homes/l/ltshu/Simulation/Rattle/struc_ase_rattle_540atoms_N200'
    out_path = os.path.join(out_dir, filename)
    write(out_path, unique_local_struc[i])