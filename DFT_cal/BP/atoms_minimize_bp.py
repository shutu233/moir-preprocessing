def Atoms_minimize(atoms):
    '''
    Args: 
    atoms: initial atoms (Atoms form)
    Returns:
    relaxed_atoms: atoms after relaxation
    eng: potential energy of atoms after relaxation(local minima)
    '''
    # Ref Johnathan's calculator example

    from ase.io import read, write  # Import read and write functions from ASE for handling atomic structures
    import numpy as np  # Import numpy for numerical operations
    from moirecompare.calculators import  NLayerCalculator  # Import custom calculators from moirecompare package
    from moirecompare.calculators import (MonolayerLammpsCalculator, 
                                        InterlayerLammpsCalculator) # Import custom calculators from moirecompare package
    from ase.optimize import FIRE, BFGS  # Import optimization algorithms from ASE

    relaxed_atoms = atoms.copy()

    # two layers
    at_1 = atoms.copy()[atoms.arrays["atom_types"] < 4]
    at_2 = atoms.copy()[atoms.arrays["atom_types"] >= 4]

    # Define atomic symbols for each layer
    layer_symbols = [["P", "P", "P", "P"],
                     ["P", "P", "P", "P"]]
    # intralayer calculator             
    intralayer_calcs = [
        MonolayerLammpsCalculator(atoms[atoms.arrays['atom_types'] < 4],
                                  layer_symbols[0],
                                  system_type='BP',
                                  intra_potential='bp.sw')
    ]

    # Initialize a list of interlayer calculators
    interlayer_calcs = []

    # Loop through the layers and set up calculators
    for i in np.arange(1, len(layer_symbols)):
        layer_atoms = atoms[
            np.logical_and(atoms.arrays['atom_types'] >= i * 4,
                           atoms.arrays['atom_types'] < (i + 1) * 4)
        ]
        # print(np.unique(layer_atoms.arrays['atom_types']))  # Print unique atom types in the current layer
        intralayer_calcs.append(MonolayerLammpsCalculator(layer_atoms,
                                                          layer_symbols=layer_symbols[i],
                                                          system_type='BP',
                                                          intra_potential='bp.sw'))

        bilayer_atoms = atoms[np.logical_and(atoms.arrays['atom_types'] >= (i - 1) * 4,
                                             atoms.arrays['atom_types'] < (i + 1) * 4)]
        # print(np.unique(bilayer_atoms.arrays['atom_types']))  # Print unique atom types in the bilayer
        # print(layer_symbols[i - 1:i + 1])  # Print symbols for the current bilayer
        interlayer_calcs.append(
            InterlayerLammpsCalculator(bilayer_atoms,
                                       layer_symbols=layer_symbols[i - 1:i + 1],
                                       system_type='BP'))
    
    # Combine the intra- and interlayer calculators into an NLayerCalculator
    n_layer_calc = NLayerCalculator(atoms,
                                    intralayer_calcs,
                                    interlayer_calcs,
                                    layer_symbols)
    atoms.calc = n_layer_calc

    # calculate the energy before relaxation
    atoms.calc.calculate(atoms)
    tot_eng = atoms.calc.results['energy']
    layer_eng1 = atoms.calc.results['layer_energy'][0][0]
    layer_eng2 = atoms.calc.results['layer_energy'][1][1]
    inter_eng = atoms.calc.results['layer_energy'][0][1]
    eng_before = [tot_eng, layer_eng1, layer_eng2, inter_eng]

    # Set up the FIRE optimizer for structural relaxation
    dyn = FIRE(atoms, trajectory=f'out.traj')
    dyn.run(fmax=1e-4)

    # Set the positions of relaxed structure (didn't change lattice constant)
    new_pos = atoms.positions
    relaxed_atoms.set_positions(new_pos)

    # calculate the energy after relaxation
    relaxed_tot_eng = atoms.calc.results['energy']
    relaxed_layer_eng1 = atoms.calc.results['layer_energy'][0][0]
    relaxed_layer_eng2 = atoms.calc.results['layer_energy'][1][1]
    relaxed_inter_eng = atoms.calc.results['layer_energy'][0][1]
    eng_after = [relaxed_tot_eng, relaxed_layer_eng1, relaxed_layer_eng2, relaxed_inter_eng]

    eng = [eng_before, eng_after]
    print(f"Relaxed energy: {atoms.calc.results['energy']}")
    return relaxed_atoms, eng


if __name__ =='__main__':
    from ase.io import read, write
    at = read('/global/homes/l/ltshu/Simulation/validate/BP/1D_BP_bilayer_twisted_relaxed_ML.xyz', format = 'extxyz')
    x = Atoms_minimize(at)
    # write("1D_BP_bilayer_twisted_relaxed.xyz", x[0], format='extxyz')

    # output the trajectory
    from ase.io.trajectory import Trajectory
    traj_path = f'out.traj'
    traj = Trajectory(traj_path)    
    images = [atom for atom in traj] 
    write(f"1D_BP_bilayer_twisted_fromML.traj.xyz", images, format="extxyz")

    # print the energy
    print(x[1])







    



