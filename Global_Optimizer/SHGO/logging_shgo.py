import logging
from logging import (
    StreamHandler,
    FileHandler,
    Formatter
)

from ase.io import read, write
from ase import Atoms
import pot_swkc
import numpy as np
from scipy.optimize import shgo
import os

logger = logging.getLogger('shgo_logger')
logger.setLevel(10)  

fh: FileHandler = FileHandler('shgo_log_540atoms.txt')
formatter: Formatter = Formatter("%(asctime)s - %(levelname)s - %(message)s")
fh.setFormatter(formatter)
fh.setLevel(10) 
logger.addHandler(fh)

xyz_file_path = "/global/homes/l/ltshu/Simulation/Structures/01-4.00_deg-540_atoms.xyz"

struc = read(xyz_file_path)

# Positions of atoms
pos = struc.positions
pos1 = np.reshape(pos, (3*len(pos), 1))

# set the bounds
# firstly set the positions \pe 0.3
posmin = pos1 - 1
posmax = pos1 + 1
bounds_array = list(zip(posmin, posmax))
bounds = [(a.item(), b.item()) for a, b in bounds_array]


def write_string(posis):
    # input: positions of atoms (format: Nx3); return: new Atoms
    struc1 = struc.copy()
    struc1.set_positions(posis)
    
    return struc1

# define the objective function
def pot_fun(arr):
    # input: atom positions (format: 3Nx1); return: potential energy
    
    arr_p = np.reshape(arr, (len(pos), 3))
    # turn atom positions array into .xyz file
    xyz = write_string(arr_p)
    pot = pot_swkc.Pot_swkc(xyz)
    return pot

# define callback function
def callback(xk):
    logger.info(f"{pot_fun(xk)}")
    #logger.info(f'Current solution: {xk}')

# shgo calculation!
res = shgo(pot_fun, bounds, n=30, callback=callback, options={"disp": True, "maxiter": 1, "seed": 100})
# res = shgo(pot_fun, bounds, n=1, callback=callback, options={"disp": True, "maxiter": 1, "maxfev": 1, "maxev": 1})

logger.info(f'{res}')
logger.info(f'{res.x}')
logger.info(f'{res.xl}')

for i in range(0, len(res.xl)):
    new_at = struc.copy()
    new_pos = np.reshape(res.xl[i], (len(pos), 3))
    new_at.set_positions(new_pos)
    filename = f'local_{i}.xyz'
    out_dir = '/global/homes/l/ltshu/Simulation/SHGO/shgo_540atoms_bo1'
    out_path = os.path.join(out_dir, filename)
    write(out_path, new_at, format='extxyz')