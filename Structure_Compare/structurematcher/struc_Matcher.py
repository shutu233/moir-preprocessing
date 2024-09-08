def structure_matcher(atoms1, atoms2):
    """
    Args:
    atoms1: 1st atoms 
    atoms2: 2nd atoms
    (in ase Atoms form)

    Returns:
    Using pymatgen.analysis.structure_matcher.StructureMatcher

    rms displacement normalized by (Vol / nsites) ** (1/3)
    and maximum distance between paired sites. 
    If no matching lattice is found None is returned.

    """
    from pymatgen.core import Structure, Lattice, Molecule
    from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.io import read, write

    struc1 = AseAtomsAdaptor.get_structure(atoms1)
    struc2 = AseAtomsAdaptor.get_structure(atoms2)

    sm = StructureMatcher(
    ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=True,
    attempt_supercell=False, comparator=ElementComparator()
    )
    # use get_rms_dist function
    rms = sm.get_rms_dist(struc1, struc2)
    
    if rms is None:
        return 100, 100
    else:
        return rms[0], rms[1]
