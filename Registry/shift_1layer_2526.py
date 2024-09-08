import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.io.ase import AseAtomsAdaptor
import os
from ase.io import read, write

# Shift rattle method
import shift_rattle

# get the relaxed structure and eng
import atoms_minimize_bei

import sys
sys.path.append('/global/homes/l/ltshu/Simulation/SHGO')
import pot_swkc
import pot_swkc_min

sys.path.append('/global/homes/l/ltshu/Simulation/struc_compare/structurematcher')
import struc_Matcher

# objective structure and unit cell of the to be moved layer
unit_cell = read('/global/homes/l/ltshu/Simulation/struc_compare/build_struc/MoS2_monolayer.xyz')
at = read('/global/homes/l/ltshu/Simulation/struc_compare/config_his/1D_25cells_0deg.xyz')

unit_vectors = unit_cell.get_cell()

# Define shifted vectors
N = 10
shifted_vectors = [] # totally 

for i in range(0, N+1):
    for j in range(0, N+1):
        shifted_vectors.append(i/N * unit_vectors[0] + j/N * unit_vectors[1])


# get the relaxed energy and relaxed shifted structure
shifted_unrelaxed_engs = []
shifted_relaxed_engs = []
shifted_relaxed_strucs = []

for i in range(0, len(shifted_vectors)):
    shifted_vector = shifted_vectors[i]
    # get the unrelaxed shifted structure
    shifted_struc = shift_rattle.Shift_rattle(at, unit_cell, shifted_vector)
    
    unrelaxed_eng = pot_swkc.Pot_swkc(shifted_struc)
    relaxed = atoms_minimize_bei.Atoms_minimize(shifted_struc)
    relaxed_eng = relaxed[1]
    relaxed_struc = relaxed[0]

    shifted_unrelaxed_engs.append(unrelaxed_eng)
    shifted_relaxed_engs.append(relaxed_eng)
    shifted_relaxed_strucs.append(relaxed_struc)

for i in range(0, len(shifted_relaxed_strucs)):
    filename = f'N10_local_{i}.xyz'
    out_dir = '/global/homes/l/ltshu/Simulation/Rattle/2526MoS21D/N10'
    out_path = os.path.join(out_dir, filename)
    write(out_path, shifted_relaxed_strucs[i], format='extxyz')



