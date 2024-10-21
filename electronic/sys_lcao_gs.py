# A SCRIPT to do  ground state  SCF calculation with LCAO basis set
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.io import vasp
import numpy as np
from ase.optimize import BFGS
from ase.io import read, write, Trajectory
from ase.parallel import paropen
system_name='ben4_ml_config2_LCAO'    # mention the system name here
struct=vasp.read_vasp("POSCAR")
#kp=[8,12,16,20,24]
#for i in range(5):
kpoint_tuple=(8,8,1)                 # k-points mention here
#write('initial.xyz',struct)
struct.pbc=[True,True,False]   # only if required! for dipolecorrection make the vaccum direction non-periodic
calc = GPAW(mode='lcao',
                 basis = 'dzp',
                 nbands='110%',
                 poissonsolver={'dipolelayer':'xy'}, # only if required! dipolecorrection in vaccum region
                 xc='vdW-DF-cx',
                 occupations=FermiDirac(0.01),
                 kpts=kpoint_tuple,
                 txt= system_name+'_gs.log',
                 parallel=dict(band=1,  # band parallelization
                               augment_grids=True,  # use all cores for XC/Poisson
                               sl_auto=True,  # enable ScaLAPACK parallelization
                               use_elpa=True))  # enable Elpa eigensolver
struct.calc = calc
energy = struct.get_potential_energy()
#ben4.calc = calc
#ben4.get_potential_energy()
ef = calc.get_fermi_level()
calc.write(system_name+'_gs.gpw')                    # if you want to continue for reuse this calculation, save the data 
#print('{0}-Total Energy = {1} '.format(i,energy))
