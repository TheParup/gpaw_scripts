# web-page: bandstructure.png
"""Band structure tutorial

Calculate the band structure of Si along high symmetry directions
Brillouin zone
"""
# P1
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.io import vasp
import numpy
from ase.build import make_supercell
#Na_supercell = make_supercell(Na_unitcell, multiplier)
# Perform standard ground state calculation (with plane wave basis)
ben4=vasp.read_vasp("POSCAR")
####
'''
#multiplier = numpy.identity(3) * 2
#ben4_sc=make_supercell(ben4,multiplier)
calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 8, 8),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01),
            txt='ben4_gs.txt')
ben4.calc = calc
ben4.get_potential_energy()
ef = calc.get_fermi_level()
calc.write('ben4_gs.gpw')
# P2
# Restart from ground state and fix potential:
calc = GPAW('ben4_gs.gpw').fixed_density(
    nbands=16,
    symmetry='off',
    kpts={'path': 'GX', 'npoints': 60},
    convergence={'bands': 8})
'''
#from gpaw import GPAW
from gpaw.unfold import Unfold, find_K_from_k

#a = 3.184
ben4_uc=vasp.read_vasp('POSCAR_uc')
PC = ben4_uc.get_cell(complete=True)
bp = PC.get_bravais_lattice().bandpath('GX', npoints=40)
x, X, _ = bp.get_linear_kpoint_axis()
M = [[4, 0, 0], [0, 4, 0], [0, 0, 1]]
Kpts = []
for k in bp.kpts:
    K = find_K_from_k(k, M)[0]
    print(K)
    Kpts.append(K)
print(len(Kpts))
# P3
#bs = calc.band_structure()
#bs.plot(filename='bandstructure.png', show=False, emax=12.0)
