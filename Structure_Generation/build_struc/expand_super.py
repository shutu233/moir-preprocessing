from ase import Atoms
from ase.build import make_supercell
from ase.io import read, write

original_structure = read('/global/homes/l/ltshu/Simulation/struc_compare/build_struc/MoS2-Bilayer_AA_6atoms.xyz')

expansion_matrix = [[10, 0, 0],
                    [0, 10, 0],
                    [0, 0, 1]]

supercell_structure = make_supercell(original_structure, expansion_matrix)

write('MoS2-Bilayer_AA_10x10.xyz', supercell_structure)
