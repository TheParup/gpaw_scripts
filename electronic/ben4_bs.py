# web-page: bandstructure.png
"""Band structure tutorial

Calculate the band structure of Si along high symmetry directions
Brillouin zone
"""
# P1
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.io import vasp
import numpy as np
from ase.build import make_supercell
#Na_supercell = make_supercell(Na_unitcell, multiplier)
# Perform standard ground state calculation (with plane wave basis)
ben4=vasp.read_vasp("POSCAR_uc")
'''
multiplier = np.identity(3) * 2
ben4_sc=make_supercell(ben4,multiplier)
#calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 7, 1),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01),
            txt='ben4_gs.txt')
#ben4.calc = calc
#ben4.get_potential_energy()
#ef = calc.get_fermi_level()
#calc.write('ben4_gs.gpw')
# P2

def read_file_to_list_of_lists(file_name):
  with open(file_name, 'r') as f:
    kptdata = []
    for line in f:
      kptdata.append(np.array(line.strip().split()[0:3],dtype=float))
  return kptdata

file_name = 'k-points'
 # print(data)
try:
  # Open the file
  f = open(file_name, 'r')
except FileNotFoundError:
  # If the file is not found, print an error message and exit the program
  print("Error: File not found")
  exit(1)

# Read the contents of the file
kptdata= read_file_to_list_of_lists(file_name)
f.close()
'''
path=ben4.cell.bandpath('GXSX1YG', npoints=200)
# -------------------------------------------------------------
# Bandstructure
# -------------------------------------------------------------
# Restart from ground state and fix potential:
calc = GPAW('ben4_gs.gpw').fixed_density(
    nbands=16,
    symmetry='off',
    kpts=path,
    convergence={'bands': 12})

# P3
bs = calc.band_structure()
bs.plot(filename='band_ORCC.png', show=False, emax=10.0)
