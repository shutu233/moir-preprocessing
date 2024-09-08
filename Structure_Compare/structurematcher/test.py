import struc_Matcher

import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.io.cif import CifWriter, CifParser
from pymatgen.io.ase import AseAtomsAdaptor

from ase.io import read, write


# for bilayer, remain layer2,translate layer1
def build_trans_structure(atoms, unitcell, N):
    """
    Args:
    atoms: the bilayer
    unitcell: monolayer, to get the Bravais vectors
    N: translate lattice constant a1, a2 / N

    return: tran_structures
    """
    
    Bra_vector = unitcell.get_cell()
    Bra_vector[2] = [0, 0, 0]
    build_trans_structure = []

    for i in range(N): # da1
        for j in range(N): # da2
            t = (i/N, j/N, 0)
            atoms2 = atoms.copy()
            atoms2.positions[atoms.arrays['atom_types'] < 3] += np.dot(t, Bra_vector)
            build_trans_structure.append(atoms2)
    return build_trans_structure

at = read("1D_25cells_0deg.xyz")
ce = read("MoS2_monolayer.xyz")

N = 2

x = build_trans_structure(at, ce, N)

write("111.xyz", x[1])
write("222.xyz", x[2])
write("333.xyz", x[3])