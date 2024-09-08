def Shift_rattle(atoms, unitcell, shifted_vector):
    '''
    Args:
    atoms: the initial bilayer structure
    unitcell: unit cell of the shifted monolayer, to get the Bravais vectors
    shifted_vector: the shifted direction & length of one target layer (shifted_vector[2] must equal to 0)
    Returns:
    shifted bilayer structure
    '''

    from ase.io import read, write
    import numpy as np

    at = atoms.copy()

    bra_vector = unitcell.get_cell()
    bra_vector[2] = [0, 0, 0]

    # projection on two Bracais vectors (as the bases)
    v1 = bra_vector[0][:2]
    v2 = bra_vector[1][:2]
    matrix = np.array([v1, v2])
    sol = np.linalg.solve(matrix.T, shifted_vector[:2])
    x = sol[0]
    y = sol[1]

    at.positions[atoms.arrays['atom_types'] < 3] += x * bra_vector[0] + y * bra_vector[1] 
    pos = at.get_positions()

    shifted_struc = atoms.copy()
    shifted_struc.set_positions(pos)

    return shifted_struc