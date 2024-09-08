def Rattle(xyz_file_path, N):
    '''
    Args:
    xyz_file_path: .xyz file of the structure to be optimized
    N: # of rattling (size of database)
    Returns: 
    all the rattled structures

    using rattle function in ase
    '''

    from ase.io import read, write
    import numpy as np
    
    at = read(xyz_file_path)

    rattled_structures = []
    for i in range(0, N):
        rattle_at = at.copy()
        rattle_at.rattle(stdev=0.2, rng=np.random)
        rattled_structures.append(rattle_at)

    return rattled_structures 
