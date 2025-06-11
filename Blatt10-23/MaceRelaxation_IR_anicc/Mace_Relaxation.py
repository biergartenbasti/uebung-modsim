# coding: utf-8
'''
script from pascal/jarno, small changes by sebastion
'''
import os
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

from ase.data.pubchem import pubchem_atoms_search #, pubchem_atoms_conformer_search # is used in last line
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.io import read, write
from ase.constraints import ExpCellFilter # only needed for periodic cells, not single molecules

from mace.calculators import mace_off
from mace.calculators import mace_mp
from mace.calculators import mace_anicc

import time

#%%
def mace_relaxation(structure, calc, fix_cell = True, f = 0.01):
    # define calculator
    structure.calc = calc
    
    # check for constraints and run
    if fix_cell == True: # kleine Anpassung für periodische Systeme, da wird ein ExpCellFilter gebraucht
        cell_filter = ExpCellFilter(structure)
        opt = BFGS(cell_filter)
    else:    
        opt = BFGS(structure)
        opt.run(fmax=f)  # Optimize until forces < 0.01 eV/Å   
        pass
    
    print("Optimizing geometry...")
    opt.run(fmax=f)  # Optimize until forces < 0.01 eV/Å
    # set up dict with relaxed structure
    structure_data={'relaxed_structure':structure,'total_energy':structure.get_total_energy()}
    return structure_data

# relaxation und schwingungen berechnen von molekül, eigentlich von Claude4.0 generiert
def calculate_ir_spectrum(structure, calc, relax=True, fix_cell = True, model_path="small"):
    """
    Calculate IR spectrum for a molecule using MACE-OFF
    
    Parameters:
    molecule: ASE Atoms object
    model_path: MACE model size ("small", "medium", "large")
    """

    # Der Teil hier kann vorerst übersprungen werden
    # Set up MACE calculator
    #calc = mace_off(model=model_path, device="cpu")  # Use "cuda" if GPU available; deaktiviert, weil ich das gerne außerhalb definieren möchte
    structure.calc = calc
    if fix_cell == True: # kleine Anpassung für periodische Systeme, da wird ein ExpCellFilter gebraucht
        cell_filter = ExpCellFilter(structure)
        pass
    
    # Optimize geometry first
    # theoretisch kann die Funktion angepasst werden und benutzt werden, Struktur wird aber über obere Funktionen relaxaiert
    if relax == True:
        print("Optimizing geometry...")
        opt = BFGS(cell_filter if fix_cell == True else structure)
        opt.run(fmax=0.01)  # Optimize until forces < 0.01 eV/Å
    else:
        print("Skipping geometry optimisation.")
        pass
    
    # Calculate vibrational frequencies
    print("Calculating vibrational frequencies...")
    vib = Vibrations(structure, name="vibrations")
    vib.run()
    
    # Get frequencies and intensities
    frequencies = vib.get_frequencies()
    
    # Print diagnostic information
    print(f"Raw frequencies shape: {frequencies.shape}")
    print(f"Frequency data types: {[type(f) for f in frequencies[:5]]}")
    
    # More robust filtering of frequencies
    real_positive_frequencies = []
    for freq in frequencies:
        if np.iscomplex(freq):
            real_part = np.real(freq)
            imag_part = np.imag(freq)
            if abs(imag_part) < 1e-6 and real_part > 0:  # Essentially real and positive
                real_positive_frequencies.append(real_part)
        elif np.isreal(freq) and freq > 0:
            real_positive_frequencies.append(float(freq))
    
    real_positive_frequencies = np.array(real_positive_frequencies, dtype=float)
    structure_data={'relaxed_structure':structure,'total_energy':structure.get_total_energy()}
    return real_positive_frequencies, vib, structure_data

# hier sollte spektrum geplottet werden, auch Claude4.0, außer savefig
def plot_ir_spectrum(frequencies, output_name='ir_spectrum', broadening=20):
    """
    Plot IR spectrum with Gaussian broadening
    
    Parameters:
    frequencies: array of vibrational frequencies in cm⁻¹
    broadening: Gaussian broadening parameter in cm⁻¹
    """
    
    # Filter out imaginary/complex frequencies and ensure real values
    real_frequencies = []
    for freq in frequencies:
        if np.iscomplex(freq):
            # Take real part if complex
            real_freq = np.real(freq)
            if real_freq > 0:  # Only include positive real frequencies
                real_frequencies.append(real_freq)
        elif np.isreal(freq) and freq > 0:
            real_frequencies.append(float(freq))
    
    real_frequencies = np.array(real_frequencies, dtype=float)
    
    print(f"Using {len(real_frequencies)} real positive frequencies out of {len(frequencies)} total")
    
    if len(real_frequencies) == 0:
        print("Warning: No positive real frequencies found!")
        return None, None
    
    # Create frequency range for plotting
    freq_range = np.linspace(500, 4000, 3500)
    spectrum = np.zeros_like(freq_range, dtype=float)
    
    # Add Gaussian peaks for each frequency
    for freq in real_frequencies:
        if 500 <= freq <= 4000:  # Typical IR range
            intensity = np.exp(-((freq_range - freq) / broadening)**2)
            spectrum += intensity
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(freq_range, spectrum, 'b-', linewidth=1.5)
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Intensity (arbitrary units)')
    plt.title(output_name+' IR Spectrum')
    plt.grid(True, alpha=0.3)
    plt.gca().invert_xaxis()  # IR spectra typically shown with decreasing wavenumber
    plt.tight_layout()
    plt.savefig(output_name+'.png', format="png")
    plt.show()
    
    return freq_range, spectrum


# Use 'multiprocessing.cpu_count()' to determine the number of available CPU cores.
cpu_count = multiprocessing.cpu_count()

# Print the number of CPU cores available on the system.
print(cpu_count)
# Step 0: Import variables from bash env, stört hier vielleicht eher
for key,value in os.environ.items():
    globals()[key]=value
    pass


#%% Relaxation and IR spectrum of Acetale input structure and Acetale Pubchem reference
 
'''
Hier ist der tatsächliche Anfang. Die Funktionen könnten perspektiFisch ausgelagert werden.
'''
t_start=time.time()

mace_calc = mace_anicc() # mace_off is an organic foundation model
molecule_name='Acetale' # wird fpr import und reference gebraucht
# Step 1: Load structure from file
# input_struc könnte auch ein input() sein, evtl einfacher
input_struc = read(molecule_name+'_in.xyz', format='xyz') # ADAPTION could be adapted to the filename specified in job script

# get reference from pubchem
reference_molecule=pubchem_atoms_search(molecule_name)
reference_molecule.calc=mace_calc
write(molecule_name+'_Pubchem.xyz',reference_molecule,format='xyz') # here die Struktur von der Pubchem database zum Vergleich

# hier werden jetzt G'schichten ausm IR-Garten geplottet
t_start=time.time()
frequencies, vib_obj, structure = calculate_ir_spectrum(input_struc, calc = mace_calc, fix_cell = False) # relaxiertes Molekül
vib_obj.clean()
t_end = time.time() - t_start 
print(f'time is {t_end} seconds.')
write('relaxed_structure.xyz', structure['relaxed_structure'],format='xyz')

frequencies_pubchem, vib_obj_pubchem, structure_reference = calculate_ir_spectrum(reference_molecule, calc = mace_calc, relax = False, fix_cell = False) # Pubchem Molekül
vib_obj_pubchem.clean()

print(f"Total free energy of relaxed structure is calculated as {structure['total_energy']} and of the Pubchem reference structure is {structure_reference['total_energy']}.") # Einheit der Energie ist bei organischen Molekülen unklar

# Print frequencies
'''
print(f"\nVibrational frequencies (cm⁻¹):")
for i, freq in enumerate(frequencies):
    print(f"Mode {i+1}: {freq:.1f} cm⁻¹")
    pass
print(f"\nVibrational frequencies (cm⁻¹):")
for i, freq in enumerate(frequencies_pubchem):
    print(f"Mode {i+1}: {freq:.1f} cm⁻¹")
    pass
'''

plotA, plotB = plot_ir_spectrum(frequencies,output_name=molecule_name, broadening=20)
plotA_pubchem, plotB_pubchem = plot_ir_spectrum(frequencies_pubchem,output_name=molecule_name+'_reference', broadening=20)
t_end = time.time() - t_start 
print(t_end)

#%% Relaxation of surface
'''
input_struc = read('POSCAR', format='vasp') # ADAPTION could be adapted to the filename specified in job script
mace_calc = mace_mp(model="small", dispersion=False, default_dtype="float32")

# relaxieren und struktur speichern
t_start=time.time()
slab = mace_relaxation(input_struc, calc = mace_calc, fix_cell = True) # relax slab
t_end = time.time() - t_start 

print(t_end)
print(f"Total free energy of relaxed structure is calculated as {slab['total_energy']}")
write('relaxed_slab.vasp', slab['relaxed_structure'],format='vasp')
'''