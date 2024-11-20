from gpaw import GPAW
import numpy as np
from ase.io import read
from collections import defaultdict, Counter
import pickle
import cProfile
def sum_pdos(doscalc_obj, energy_arr, atom_idx_list, orbital_idx, m_idx=None):
    # Initialize the list to store PDOS for each atom
    pdos_orbs = []

    # Default handling for no magnetic quantum number provided
    if m_idx is None:
        for i in atom_idx_list:
            try:
                element_orb = doscalc_obj.raw_pdos(energy_arr, a=i, l=orbital_idx, width=0.0)
                pdos_orbs.append(element_orb)
            except Exception as e:
                print(f"Error processing atom index {i} without m_idx: {e}")
    # Handling for specific magnetic quantum number provided
    else:
        for i in atom_idx_list:
            try:
                element_orb = doscalc_obj.raw_pdos(energy_arr, a=i, l=orbital_idx, m=m_idx, width=0.0)
                pdos_orbs.append(element_orb)
            except Exception as e:
                print(f"Error processing atom index {i} with m_idx: {e}")

    # Check if pdos_orbs is not empty before summing
    if not pdos_orbs:
        raise ValueError("The pdos_orbs list is empty. Ensure that the atom indices and orbital indices are correct.")

    # Perform the summation along the first axis
    pdos = np.sum(pdos_orbs, axis=0)
    return pdos

def atoms_metadata(atom_obj):
    # Initialize a dictionary for metadata
    metadata = {}

    # Dictionary to keep track of counts and indices of each species
    species_info = defaultdict(lambda: {'count': 0, 'indices': []})
    # Populate species_info dictionary
    for idx, atom in enumerate(atom_obj):
        symbol = getattr(atom, 'symbol', None)
        if symbol is None:
            print(f"Skipping atom at index {idx} without symbol attribute.")
            continue
        species_info[symbol]['count'] += 1
        species_info[symbol]['indices'].append(idx)

    # Populate metadata dictionary
    atoms_list = atom_obj.get_chemical_symbols()
    species_list = list(Counter(atoms_list))
    metadata['unique_species'] = species_list
    metadata['num_unique_species'] = len(species_info)
    metadata['species_info'] = dict(species_info)  # Convert back to regular dictionary for better readability

    return metadata

####----MAIN PROGRAM-----###########################
def main():
    name = 'ben4_NR_5x5_N_vac_clean'  # Change input name here
    filename = f'{name}_gs.gpw'
    emin = -10
    emax = 10
    res = 0.01

    # Calculate the number of energy bins
    bins = int((emax - emin) / res) + 1

    # Load the structure and calculation
    calc_struct = GPAW(filename, txt=None)
    struct = read(filename)
    # Get metadata about the atomic structure
    struct_data = atoms_metadata(struct)

    # Calculate the density of states
    struct_doscalc = calc_struct.dos()
    energy_arr = np.linspace(emin, emax, bins)
    total_dos = struct_doscalc.raw_dos(energy_arr,width=0.0)
    # Process each unique species in the structure
    species_list = struct_data['unique_species']
    # Dictionary to keep track of (shell/oribital) projected dos data for each atomic species
    pdos_orbital_data = {} #defaultdict(lambda: {'s': np.zeros(energy_arr.size), 'p': np.zeros(energy_arr.size)})
    #Loop over each species for pdos of species' orbitals
    for atm in species_list:
        atm_idx_list = struct_data['species_info'][atm]['indices']
        # Compute PDOS for s and p orbitals
        pdos_s = sum_pdos(struct_doscalc, energy_arr, atm_idx_list, orbital_idx=0)
        pdos_px = sum_pdos(struct_doscalc, energy_arr, atm_idx_list, orbital_idx=1,m_idx=2) #look documentation for m index
        pdos_py = sum_pdos(struct_doscalc, energy_arr, atm_idx_list, orbital_idx=1,m_idx=0)
        pdos_pz = sum_pdos(struct_doscalc, energy_arr, atm_idx_list, orbital_idx=1,m_idx=1)
        #initize a dictionary and dump pdos for each unique element use pickle to dump
        temp_dict={'s':pdos_s,'px':pdos_px,'py':pdos_py,'pz':pdos_pz}
        #
        pdos_orbital_data[atm]=temp_dict
        #pdos_shell_data[atm]['p']=pdos_p
        # species_info[atm]['p'].append(idx)
    #----END-of---LOOP -----### 
    dos_data = {'energies':energy_arr,
                 'dos_data':total_dos,
                 'pdos_data':pdos_orbital_data }
    # Save the band structure data to a pickle file

    with open('{0}_pdos_data.pkl'.format(name), 'wb') as f:

              pickle.dump(dos_data, f)
    '''
        # Initialize and fill the PDOS data array
        pdos_data = np.zeros((energy_arr.size, 3))
        pdos_data[:, 0] = energy_arr
        pdos_data[:, 1] = pdos_s
        pdos_data[:, 2] = pdos_p

        # Define column names and save data to a text file
        column_names = f"energy {atm}-s {atm}-p"
        np.savetxt(f'{atm}_shell_pdos.dat', pdos_data, delimiter=' ', header=column_names, comments='')
      '''
#   print('Done')
if __name__ == "__main__":
    main()
    #cProfile.run('main()') # if you want run with cProfiling uncomment here and comment previous line i.e main()
