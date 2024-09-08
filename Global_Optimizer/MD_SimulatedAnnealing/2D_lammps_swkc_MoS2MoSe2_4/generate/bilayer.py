#!/usr/bin/env python

# Generates a set of bilayers MoS2/MoSe2 with varying twist angles from 0 to
# 60°.
#
# Felipe H. da Jornada (2020)


import numpy as np
import pymoire as pm

p = pm.materials.get_materials_db_path()
#p = "/global/homes/l/ltshu/python/pymoire/pymoire/pymoire/c2db_structures"

layer1 = pm.read_monolayer(p/'MoS2.cif')
layer2 = pm.read_monolayer(p/'MoSe2.cif')

# Attempt to generate twisted structures with all possible angles, from 0° to
# 60° (including), and with a spacing of 2°. The dθ/10 term ensures that the
# end point is included.
dθ = .2
# angles = np.arange(0, 60 + dθ/10, dθ)
angles = np.arange(0, 5 + dθ/10, dθ)
# angles = np.array([20.,])

pm.gen_moire(layer1, layer2, angles) #, tol_strain=0.1, R_max_xtal=2)
