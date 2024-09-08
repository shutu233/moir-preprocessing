# output the trajectory
from ase.io.trajectory import Trajectory
from ase.io import write
traj_path = f'out.traj'
traj = Trajectory(traj_path)    
images = [atom for atom in traj] 
write(f"out.traj.xyz", images, format="extxyz")