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
                                        InterlayerLammpsCalculator) # Import custom calculators from moirecompare package
    from ase.optimize import FIRE, BFGS  # Import optimization algorithms from ASE

    relaxed_atoms = atoms.copy()

    # two layers
    at_1 = atoms.copy()[atoms.arrays["atom_types"] < 3]
    at_2 = atoms.copy()[atoms.arrays["atom_types"] >= 3]

    # Define atomic symbols for each layer
    layer_symbols = [["Mo", "S", "S"],
                     ["Mo", "Se", "Se"]]

    # intralayer calculator             
    intralayer_calcs = [
        MonolayerLammpsCalculator(atoms[atoms.arrays['atom_types'] < 3],
                                  layer_symbols[0],
                                  system_type='TMD',
                                  intra_potential='tmd.sw')
    ]

    # Initialize a list of interlayer calculators
    interlayer_calcs = []

    # Loop through the layers and set up calculators
    for i in np.arange(1, len(layer_symbols)):
        layer_atoms = atoms[
            np.logical_and(atoms.arrays['atom_types'] >= i * 3,
                           atoms.arrays['atom_types'] < (i + 1) * 3)
        ]
        # print(np.unique(layer_atoms.arrays['atom_types']))  # Print unique atom types in the current layer
        intralayer_calcs.append(MonolayerLammpsCalculator(layer_atoms,
                                                          layer_symbols=layer_symbols[i],
                                                          system_type='TMD',
                                                          intra_potential='tmd.sw'))

        bilayer_atoms = atoms[np.logical_and(atoms.arrays['atom_types'] >= (i - 1) * 3,
                                             atoms.arrays['atom_types'] < (i + 1) * 3)]
        # print(np.unique(bilayer_atoms.arrays['atom_types']))  # Print unique atom types in the bilayer
        # print(layer_symbols[i - 1:i + 1])  # Print symbols for the current bilayer
        interlayer_calcs.append(
            InterlayerLammpsCalculator(bilayer_atoms,
                                       layer_symbols=layer_symbols[i - 1:i + 1],
                                       system_type='TMD'))
    
    # Combine the intra- and interlayer calculators into an NLayerCalculator
    n_layer_calc = NLayerCalculator(atoms,
                                    intralayer_calcs,
                                    interlayer_calcs,
                                    layer_symbols)
    atoms.calc = n_layer_calc

    # calculate the energy before relaxation
    atoms.calc.calculate(atoms)
    layer_eng1 = atoms.calc.results['layer_energy'][0][0]
    layer_eng2 = atoms.calc.results['layer_energy'][1][1]
    inter_eng = atoms.calc.results['layer_energy'][0][1]
    eng_before = [layer_eng1, layer_eng2, inter_eng]
    # Set up the FIRE optimizer for structural relaxation
    dyn = FIRE(atoms)
    dyn.run(fmax=1e-4)
    # calculate the energy after relaxation
    relaxed_layer_eng1 = atoms.calc.results['layer_energy'][0][0]
    relaxed_layer_eng2 = atoms.calc.results['layer_energy'][1][1]
    relaxed_inter_eng = atoms.calc.results['layer_energy'][0][1]
    eng_after = [relaxed_layer_eng1, relaxed_layer_eng2, relaxed_inter_eng]
    return eng_before, eng_after


