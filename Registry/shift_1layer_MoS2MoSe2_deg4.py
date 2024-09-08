import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.io.ase import AseAtomsAdaptor
import os
from ase.io import read, write
from ase import Atoms

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
at = read('/global/homes/l/ltshu/Simulation/Structures/01-4.00_deg-540_atoms.xyz')
Mo1 = np.array([0.03177610, 1.79521115, 0])
Mo2 = np.array([1.67757683, 4.46886443, 0])
Mo3 = np.array([-1.51868681, 4.50772995, 0])
v1 = Mo2 - Mo1
v2 = Mo3 - Mo1

cell = [list(v1), list(v2), [0, 0, 12]]
positions = [[0.03177610, 1.79521115, 0], [1.61402464, 0.87844212, -1.58592600], [1.61402464, 0.87844212, 1.58592600]]
unit_cell = Atoms('MoS2', positions=positions, cell=cell)
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
    out_dir = '/global/homes/l/ltshu/Simulation/Rattle/540atoms_4deg/N10'
    out_path = os.path.join(out_dir, filename)
    write(out_path, shifted_relaxed_strucs[i], format='extxyz')



