# Define a function, calculate the potental energy using SW.KC potential
# Use Johnathan's calculator

def Pot_swkc(at):
    '''
    input: Atoms
    return: relaxed potential energy
    '''
    from ase.io import read, write
    from pathlib import Path
    from moirecompare.calculators import MonolayerLammpsCalculator, BilayerLammpsCalculator
    from ase.calculators.lammpslib import LAMMPSlib
    from moirecompare.calculators import QECalculator
    from ase.optimize import FIRE

    # at = read(xyz_file_path,format='extxyz')
    chem_syms = [["Mo", "S", "S"], ["Mo", "S", "S"]]

    # Calculate energy (after relaxation)
    calc = BilayerLammpsCalculator(at, chemical_symbols = chem_syms, system_type = 'TMD')
    
    at.calc = calc
    # get the energy
    # at.calc.calculate(at)
    # do relaxation
    at.calc.minimize(at)

    return at.calc.results['energy']
    




    