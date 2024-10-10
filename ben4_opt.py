#python script to ionic relax a structure using LCAO basis
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.io import vasp
import numpy as np
from ase.optimize import BFGS
from ase.io import read, write, Trajectory
from ase.parallel import paropen
system_name='ben4'
struct=vasp.read_vasp("POSCAR_new")
kpoint_tuple=(100,1,1)
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
write('final.cif',struct)
vasp.write_vasp('CONTCAR',struct)

