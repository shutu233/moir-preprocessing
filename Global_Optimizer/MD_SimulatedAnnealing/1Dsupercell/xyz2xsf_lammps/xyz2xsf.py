from ase import io
from ase.io.lammpsdata import write_lammps_data

xyzfile = io.read("MoS2WSe2_1D_Stretch_26atoms.xyz")

io.write("MoS2WSe2_1D_Stretch_26atoms.xsf", xyzfile)

# generate lammps data file, but格式和后续计算sw/kc力需要的原子种类不一定一致
write_lammps_data("MoS2WSe2_1D_Stretch_26atoms.dat", xyzfile)