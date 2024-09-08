from ase.io import read, write
import numpy as np

at = read('/global/homes/l/ltshu/Simulation/Rattle/2526MoS21D/N10/N10_local_0.xyz')
at_new = at.copy()
cell = at.get_cell()
pos = at.get_positions()
pos_new = pos.copy()

zmin = np.max(pos[:][2])
zmax = np.min(pos[:][2])
zmiddle = (zmin + zmax)/2

for i in range(0, len(pos)):
    pos_new[i][2] += cell[2][2]/2 + zmiddle

at_new.set_positions(pos_new)
write('2526_local0.xyz', at_new, format='xyz')