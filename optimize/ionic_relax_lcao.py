#python script to ionic relax a structure using LCAO basis
# INPUT File should be in vasp POSCAR FORMAT
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.io import vasp
import numpy as np
from ase.optimize import BFGS
from ase.io import read, write, Trajectory
from ase.parallel import paropen
system_name=' '                 #  MENTION system name
struct=vasp.read_vasp("POSCAR") # MENTION POSCAR FILE NAME
kpoint_tuple=(100,1,1)          # CHANGE THE KPOINTS as per requirement
write('initial.cif',struct)
calc = GPAW(mode='lcao',
            basis='dzp',
            nbands='110%',
            xc='PBE',
            occupations=FermiDirac(0.01),
            kpts=kpoint_tuple,
            parallel=dict(band=1,  # band parallelization
                          augment_grids=True,  # use all cores for XC/Poisson
                          sl_auto=True,  # enable ScaLAPACK parallelization
                          use_elpa=True))  # enable Elpa eigensolver
struct.calc = calc
#with paropen('opt.traj', "a") as fd:
opt = BFGS(struct)#, trajectory=fd)
traj = Trajectory(system_name + '_opt.traj', 'w',struct)
opt.attach(traj)
opt.run(fmax=0.03)
write('final.cif',struct)  #OPTIMIZED STRUCTURE IN CIF FILE
vasp.write_vasp('CONTCAR',struct)  # OPTIMIZED STRUCTURE IN CONTCAR FORMAT

