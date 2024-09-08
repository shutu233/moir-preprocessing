from ase import io
from ase.io.lammpsdata import write_lammps_data

xyzfile = io.read("/global/homes/l/ltshu/Simulation/moiretest/6atoms/xyz2xsf_lammps/MoS2-Bilayer_AA_6atoms.xyz")

io.write("MoS2-Bilayer_AA_6atoms.xsf", xyzfile)

# generate lammps data file, but格式和后续计算sw/kc力需要的原子种类不一定一致
write_lammps_data("MoS2-Bilayer_AA_6atoms.dat", xyzfile)