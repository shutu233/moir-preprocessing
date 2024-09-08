#!/usr/bin/env python

# Converts an xsf (or similar crystal file) to a LAMMPS input.
# Only supports TMD bilayers now with SW + KC interactions.
#
# Jornada group (2021)

import numpy as np
import ase
import ase.io
import ase.data
import sys


def unit_range_fixed(x, eps=1e-6):
    y = x.copy()
    y[(np.fabs(y)>0)&(np.fabs(y)<eps)] = 0
    y = y%1
    y[y>1-eps] = 0
    return y


def get_X(symbs):
    if 'S' in symbs:
        return 'S'
    elif 'Se' in symbs:
        return 'Se'
    elif 'Te' in symbs:
        return 'Te'
    else:
        raise ValueError(symbs)


def get_M(symbs):
    if 'Mo' in symbs:
        return 'Mo'
    elif 'W' in symbs:
        return 'W'
    else:
        raise ValueError(symbs)


fname_in = sys.argv[1]
atoms = ase.io.read(fname_in)

# need to change according to the position of two layers
cond_b = atoms.positions[:,2] < 5
cond_t = atoms.positions[:,2] >= 5

symbs_b = list(set(atoms.symbols[cond_b]))
symbs_t = list(set(atoms.symbols[cond_t]))

M_b = get_M(symbs_b)
X_b = get_X(symbs_b)
M_t = get_M(symbs_t)
X_t = get_X(symbs_t)

symbs_layer = [[M_b, X_b, X_b], [M_t, X_t, X_t]]
symbs_lammps = [M_b, X_b, X_b, M_t, X_t, X_t]
print(symbs_lammps)

avgz_b = atoms.positions[cond_b][:,2].mean()
avgz_t = atoms.positions[cond_t][:,2].mean()
#print(avgz_b)
#print(avgz_t)

tol = 1e-1
lammps_atoms = np.zeros((len(atoms),), dtype=int)
posz = atoms.positions[:,2]

cond = cond_b & (posz < avgz_b)
lammps_atoms[cond] = 2
cond = cond_b & (posz > avgz_b)
lammps_atoms[cond] = 3
cond = cond_b & (np.fabs(posz - avgz_b) < tol)
lammps_atoms[cond] = 1

cond = cond_t & (posz < avgz_t)
lammps_atoms[cond] = 5
cond = cond_t & (posz > avgz_t)
lammps_atoms[cond] = 6
cond = cond_t & (np.fabs(posz - avgz_t) < tol)
lammps_atoms[cond] = 4

#print(lammps_atoms)
assert np.all(lammps_atoms!=0)

masses = [ase.data.atomic_masses_iupac2016[
            ase.data.atomic_numbers[symb]] for symb in symbs_lammps]

pos = atoms.positions.copy()
#pos[:,0] -= np.amin(pos[:,0])
#pos[:,1] -= np.amin(pos[:,1])
pos[:,2] -= np.amin(pos[:,2])

# https://lammps.sandia.gov/doc/Howto_triclinic.html

# Rewrite positions in crystal coords
# pos = pos_crys * cell
pos_crys = pos @ np.linalg.inv(atoms.cell)
pos_cys = unit_range_fixed(pos_crys)
print(pos_crys)

# LAMMPS requires the cell in a specific way. In Python order:
# [ax,  0,  0]
# [bx, by,  0]
# [cx, cy, cz]
cell = atoms.cell.copy()
print(cell)
Q, cell2 = np.linalg.qr(cell.T)
cell2 = cell2.T

print(cell2)
#print(cell2)
if cell2[0,0] < 0:
    cell2[:,:2] *= -1
print(cell2)
if cell2[1,1] < 0:
    cell2[:,1] *= -1
    cell2[:,2] *= -1
print(cell2)
if cell2[2,2] < 0:
    cell2[2,2] *= -1
#assert np.linalg.det(cell2) > 0
print(cell2)

assert np.allclose(
        np.linalg.norm(cell, axis=1),
        np.linalg.norm(cell2, axis=1))

pos = pos_cys @ cell2
#pos = pos @ Q
#print(pos)

# delete the periodic structure
atoms2 = atoms.copy()
atoms2.cell = cell2
atoms2.positions = pos
# z direction
atoms2.positions[:,2] += 2
atoms2.write('lammps_in.xsf')

with open('lammps.dat', 'w') as f:
    f.write(f'''
 {len(atoms)} atoms
 6 atom types

 0.000000 {cell2[0,0]} xlo xhi
 0.000000 {cell2[1,1]} ylo yhi
 0.000000 100.00000 zlo zhi
 {cell2[1,0]} 0.000000  0.000000 xy xz yz

 Masses 

''')
    #1 183.840000 
    data = np.column_stack((np.arange(len(masses))+1, masses))
    np.savetxt(f, data, fmt='%d %12.6f')

    f.write('''
 Atoms

''')
    #1 2 33.725809 61.747831 4.344402

    data = np.column_stack((
        np.arange(len(atoms))+1, lammps_atoms, pos))
    np.savetxt(f, data, fmt='%4d %1d %12.6f %12.6f %12.6f')


def get_SW_str(ilayer):
    l = ['NULL'] * 6
    off = (ilayer - 1) * 3
    l[off + 0] = symbs_lammps[off + 0]
    l[off + 1] = symbs_lammps[off + 1]
    l[off + 2] = symbs_lammps[off + 2]
    print(l)
    return ' '.join(l)


def get_KC_str(i1, i2):
    l = ['NULL'] * 6
    l[i1-1] = symbs_lammps[i1-1]
    l[i2-1] = symbs_lammps[i2-1]
    return ' '.join(l)


with open('lammps.in', 'w') as f:
    f.write(f'''\
#Initialize--
#general settings
units           metal
dimension       3
box tilt        large
atom_style      atomic

# structure
boundary        p p p
read_data       lammps.dat

# potentials
pair_style hybrid/overlay sw/mod sw/mod kolmogorov/crespi/z 14.0 kolmogorov/crespi/z 14.0 kolmogorov/crespi/z 14.0 kolmogorov/crespi/z 14.0 lj/cut 10.0

# Intralayer Interaction
pair_coeff * * sw/mod 1 tmd.sw {get_SW_str(1)}
pair_coeff * * sw/mod 2 tmd.sw {get_SW_str(2)}

# Interlayer Interaction
pair_coeff 1 5 kolmogorov/crespi/z 1 WS.KC {get_KC_str(1,5)}
pair_coeff 3 4 kolmogorov/crespi/z 2 WS.KC {get_KC_str(3,4)}
pair_coeff 3 5 kolmogorov/crespi/z 3 WS.KC {get_KC_str(3,5)}
pair_coeff 1 4 kolmogorov/crespi/z 4 WS.KC {get_KC_str(1,4)}
pair_coeff * * lj/cut 0.0 3.0
neighbor        2.0 bin
neigh_modify every 1 delay 0 check yes

#optimize at 0 K
dump            1 all custom 100 dump.initial id type x y z
thermo          1000
thermo_style    custom step pe press
undump          1

min_style       fire
minimize        0.0 1.0e-4 1000000 1000000
write_data      lammps.dat_min
''')