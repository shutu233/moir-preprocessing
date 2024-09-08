from ase.io import read, write
from ase import Atoms
import numpy as np

at = read('MoS2_monolayer.xyz') 

pos = at.positions
new_pos = pos + np.array([0, 0, 3.08592600]) - np.array([0, 0, -3.08592600])

new_at = at.copy()
new_at.set_positions(new_pos)

combined_at = at + new_at
post = combined_at.get_positions()
post[:, 2] += 2.0

combined_at.set_positions(post)
write('output.xyz', combined_at)
