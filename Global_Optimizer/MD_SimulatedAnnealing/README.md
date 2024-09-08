### Basic Steps of Simulated Annealing
1. Use Pymoiré/1D-gen to generate a structure
2. Turn .xyz or .xsf file into lammps data format
3. `in.xx`: input file for LAMMPS.`in.relax`: force field relaxation; `in.sa_linear`: linear cooling; `in.sa_step`: cooling step by step (all in SW+KC force field)
