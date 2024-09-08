#!/usr/bin/env python

from pathlib import Path
import pickle
import numpy as np
import matplotlib.pyplot as plt


def analyze_results(fname_config, plot=False):

    print(f'Config file: {fname_config}')
    with open(fname_config, 'rb') as f:
        dct = pickle.load(f)

    output_dir = dct['output_dir']
    output_dir = Path(fname_config).parent / Path(output_dir)
    print(f'Output dir: {output_dir}')

    plist = sorted(output_dir.glob('*.log'))
    print(f'Found {len(plist)} structures.')

    natoms = []
    strains = []
    angles = []
    fnames = []

    for p in plist: 
        with open(p.with_suffix('.pickle'), 'rb') as f:
            ret = pickle.load(f)
            strain = ret.strains.max()
            natoms.append(len(ret.structure))
            strains.append(strain)
            angles.append(ret.theta_deg)
            fnames.append(p.name)

    strains = np.array(strains)
    natoms = np.array(natoms)

    if plot:
        plt.figure(figsize=(12,6))

        colors = ['r', 'b']

        #plt.subplot(121)

        ax1 = plt.gca()
        color = colors[0]
        ax1.set_xlabel('Twist angle (°)')
        ax1.set_ylabel('Strain (%)', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.plot(angles, strains*100, 'o', ms=5, color=color)

        ax2 = ax1.twinx()
        color = colors[1]
        ax2.set_ylabel('Number of atoms', color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.plot(angles, natoms, '^', ms=5, color=color)

        #plt.subplot(122)
        #plt.plot(strains*100, natoms, '*')
        #plt.xlabel('Strain (%)')
        #plt.ylabel('Number of atoms')

        plt.tight_layout()

        plt.savefig(f'strain_vs_atoms.pdf', dpi=200, bbox_inches='tight', pad_inches=0.05)
        plt.show()


    cost = strains.copy()
    #cost[natoms > 150] = np.inf
    cost += natoms * 1e3
    idx_best = np.argmin(cost)
    print()
    print('Optimal structure:')
    print(f'- Index: {idx_best}')
    print(f'- File name: {fnames[idx_best]}')
    print(f'- Strain: {strains[idx_best]*100:.3f}%')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fname_config', default='pymoire.pickle', nargs='?')
    parser.add_argument('--plot', default=False, action='store_true')
    args = parser.parse_args()
    analyze_results(**vars(args))