def Intra_interlayer_Eng(atoms):
    '''
    Args:
    atoms: initial atoms (Atoms form, bilayer structure)
    
    Returns:

    layer_eng1: the intralayer energy of the 1st layer before relaxation
    layer_eng2: the intralayer energy of the 2nd layer before elaxation
    inter_eng: the interlayer energy between two layers before relaxation
    relaxed_xxx: the energy after relaxation

    eng_before, eng_after: (relaxed)[layer_eng1, layer_eng2, inter_eng]
    '''

    from ase.io import read, write  # Import read and write functions from ASE for handling atomic structures
    import numpy as np  # Import numpy for numerical operations
    from moirecompare.calculators import  NLayerCalculator  # Import custom calculators from moirecompare package
    from moirecompare.calculators import (MonolayerLammpsCalculator, 
                                        InterlayerLammpsCalculator, 
                                        BilayerLammpsCalculator) # Import custom calculators from moirecompare package
    from ase.optimize import FIRE, BFGS  # Import optimization algorithms from ASE

    relaxed_atoms = atoms.copy()

    # Define atomic symbols for each layer
    layer_symbols = [["Mo", "S", "S"],
                     ["Mo", "S", "S"]]

    
    # Combine the intra- and interlayer calculators into an NLayerCalculator
    bilayer_calc = BilayerLammpsCalculator(atoms, chemical_symbols=layer_symbols, system_type='TMD')
    atoms.calc = bilayer_calc

    # calculate the energy before relaxation
    atoms.calc.calculate(atoms)
    layer_eng1 = atoms.calc.results['L1_energy']
    layer_eng2 = atoms.calc.results['L2_energy']
    inter_eng = atoms.calc.results['IL_energy']
    eng_before = [layer_eng1, layer_eng2, inter_eng]
    # Set up the FIRE optimizer for structural relaxation
    dyn = FIRE(atoms, trajectory=f'out2.traj')
    dyn.run(fmax=1e-3)

    relaxed_layer_eng1 = atoms.calc.results['L1_energy']
    relaxed_layer_eng2 = atoms.calc.results['L2_energy']
    relaxed_inter_eng = atoms.calc.results['IL_energy']
    eng_after = [relaxed_layer_eng1, relaxed_layer_eng2, relaxed_inter_eng]

    # output the trajectory
    from ase.io.trajectory import Trajectory
    traj_path = f'out2.traj'
    traj = Trajectory(traj_path)    
    images = [atom for atom in traj] 
    write(f"out2.traj.xyz", images, format="extxyz")

    return eng_before, eng_after


