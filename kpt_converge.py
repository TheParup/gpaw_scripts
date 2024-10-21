# script to converge kpt interatively over loop of kpoints
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.io import vasp
import numpy as np
from ase.optimize import BFGS
from ase.io import read, write, Trajectory
from ase.parallel import paropen, parprint
system_name='ben4_clean_44'
struct=vasp.read_vasp("POSCAR")
kp=[4,6,8,10,12,16]
for i in range(6):
    kpoint_tuple=(kp[i],kp[i]-1,1)
    #write('initial.xyz',struct)
    struct.pbc=[True,True,False]
    calc = GPAW(mode='lcao',
                 basis = 'dzp',
                 nbands='110%',
                 poissonsolver={'dipolelayer':'xy'},
                 xc='vdW-DF-cx',
                 occupations=FermiDirac(0.01),
                 kpts=kpoint_tuple,
                 txt = system_name+str(kp[i])+'.log'
                 parallel=dict(band=1,  # band parallelization
                               augment_grids=True,  # use all cores for XC/Poisson
                               sl_auto=True,  # enable ScaLAPACK parallelization
                               use_elpa=True))  # enable Elpa eigensolver
    struct.calc = calc
    energy = struct.get_potential_energy()
    #ben4.calc = calc
    #ben4.get_potential_energy()
    ef = calc.get_fermi_level()
    #calc.write(system_name+'_gs.gpw')
    parprint('{0}-Total Energy = {1} '.format(i,energy))
~                                                      
