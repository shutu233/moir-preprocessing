from ase.io import read, write
at = read('/global/homes/l/ltshu/Simulation/validate/BP/bp_bp_traj.xyz', index = -1, format = 'extxyz')
at_old = read('/global/homes/l/ltshu/Simulation/validate/BP/1D_BP_bilayer_twisted.xyz', format='extxyz')

pos_new = at.get_positions()
at_new = at_old.copy()
at_new.set_positions(pos_new)
write("1D_BP_bilayer_twisted_relaxed_ML.xyz", at_new, format='extxyz')