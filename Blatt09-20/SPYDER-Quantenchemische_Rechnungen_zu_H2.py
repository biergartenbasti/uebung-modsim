#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 14:50:42 2025

@author: sutz
"""
'''
# installing modules
%pip install torch
%pip install mace-torch
%pip install ase
%pip install nglview
'''

#%% import modules

# MACE um die Rechner und Modelle zu definieren
from mace.calculators import mace_anicc # coupled cluster model
from mace.calculators import mace_mp # coupled cluster model
from mace.calculators import MACECalculator # Calculator um Rechner zu definieren

# ASE um Strukturen als atoms object zu lesen, schreiben und ändern
from ase.io import read, write # Zum Lesen und Schreiben von Strukturdateien
from ase.visualize import view # Zur Visualisierung von Strukturen
from ase.constraints import ExpCellFilter # Wird für die Strukturoptimierung von periodischen Systemen benötigt
from ase.optimize import BFGS # Solver für die Berechnungen

# nur für Utility
import copy # Objekt muss kopiert werden, sonst wird es verändert
import time # Um Zeit für den Funktionsdurchlauf zu berechnen
import numpy as np # Utility für Vektoren

#%% define relaxation
def relaxation_mace(atoms, calc, fix_cell = True, f = 0.01):
    '''
    Parameters
    ----------
    structure : ATOMS object
        Struktur die für die Optimierung genutzt werden soll.
    calc : CALCULATOR object
        Rechner der im weiteren Laufe zur Simulation genutzt werden soll.
    fix_cell : BOOL, optional
        Beschreibt, ob die Gitterkonstanten fix oder veränderlich sind. The default is True.
    f : FLOAT, optional
        Konvergenzkriterium zum Beenden der Optimierung. The default is 0.01.

    Returns
    -------
    structure : ATOMS object
        Beinhaltet die optimierte Struktur als ATOMS object.
    '''
    # define calculator
    t_start = time.time()
    structure = copy.deepcopy(atoms)
    structure.calc = calc # Hier wird ein zusätzlicher Rechner definiert, auf Grundlage dessen alle weiteren Rechnungen durchgeführt werden
    
    # check for constraints and run
    if fix_cell == True: # kleine Anpassung für periodische Systeme, da wird ein ExpCellFilter gebraucht
        cell_filter = ExpCellFilter(structure)
        opt = BFGS(cell_filter) # BFGS ist ein Solver
    else:    
        opt = BFGS(structure)
        opt.run(fmax=f)  # Optimize until forces < 0.01 eV/Å   
        pass
    
    print("Optimizing geometry...")
    opt.run(fmax=f)  # Optimise until forces < 0.01 eV/Å
    # set up dict with relaxed structure
    t_end = time.time()
    print(f'\nDie Funktion brauchte {t_end-t_start} Sekunden.')
    return structure

#%% Start calculations

if __name__ == '__main__':
    # Definiere calculators
    calc_gga = MACECalculator(model_paths="./MACE-matpes-pbe-omat-ft.model",device='cpu')
    calc_metagga = MACECalculator(model_paths="./MACE-matpes-r2scan-omat-ft.model",device='cpu')
    calc_hybrid = MACECalculator(model_paths='./MACE-OFF24_medium.model',device='cpu') # medium heißt hier, das Modell ist kleiner als der Standard
    calc_cc = mace_anicc(device='cpu')
    
    # import structures
    H2 = read('./H2.xyz',format='xyz')
    H = read('./H.xyz',format='xyz')
    #view(H2,viewer='ngl')
    
    # Simulationsteil
    # coupled cluster
    h2_cc = relaxation_mace(H2, calc_cc, fix_cell = False, f = 0.01)
    h_cc = relaxation_mace(H, calc_cc, fix_cell = False, f = 0.01)
    
    # hybrid
    h2_hybrid = relaxation_mace(H2, calc_hybrid, fix_cell = False, f = 0.01)
    h_hybrid = relaxation_mace(H, calc_hybrid, fix_cell = False, f = 0.01)
    
    # metagga
    h2_metagga = relaxation_mace(H2, calc_metagga, fix_cell = False, f = 0.01)
    h_metagga = relaxation_mace(H, calc_metagga, fix_cell = False, f = 0.01)
    
    # gga
    h2_gga = relaxation_mace(H2, calc_gga, fix_cell = False, f = 0.01)
    h_gga = relaxation_mace(H, calc_gga, fix_cell = False, f = 0.01)
    
    # print alle Atomabstände
    print(f"Der Atomabstand im Coupled Cluster Modell ist {np.linalg.norm(h2_cc.get_positions()[1]-h2_cc.get_positions()[0]):.2f} Angstrom.\n")
    print(f"Der Atomabstand im Hybrid Modell ist {np.linalg.norm(h2_hybrid.get_positions()[1]-h2_hybrid.get_positions()[0]):.2f} Angstrom.\n")
    print(f"Der Atomabstand im METAGGA Modell ist {np.linalg.norm(h2_metagga.get_positions()[1]-h2_metagga.get_positions()[0]):.2f} Angstrom.\n")
    print(f"Der Atomabstand im GGA Modell ist {np.linalg.norm(h2_gga.get_positions()[1]-h2_gga.get_positions()[0]):.2f} Angstrom.\n")
    
    view(h2_hybrid,viewer='ngl')
    
    # Bindungsenergie
    # Coupled Cluster
    Ebind_cc = h2_cc.get_total_energy() - 2*h_cc.get_total_energy()
    print(f" Die Bindungsenergie im Coupled Cluster Modell ist {Ebind_cc:.2f} eV oder {Ebind_cc*96.485:2f} kJ/mol.")
    # Hybrid
    Ebind_hybrid = h2_hybrid.get_total_energy() - 2*h_hybrid.get_total_energy()
    print(f" Die Bindungsenergie im Hybrid Modell ist {Ebind_hybrid:.2f} eV oder {Ebind_hybrid*96.485:2f} kJ/mol.")
    
    # METAGGA
    Ebind_metagga = h2_metagga.get_total_energy() - 2*h_metagga.get_total_energy()
    print(f" Die Bindungsenergie im METAGGA Modell ist {Ebind_metagga:.2f} eV oder {Ebind_metagga*96.485:2f} kJ/mol.")
    
    # GGA
    Ebind_gga = h2_gga.get_total_energy() - 2*h_gga.get_total_energy()
    print(f" Die Bindungsenergie im GGA Modell ist {Ebind_gga:.2f} eV oder {Ebind_gga*96.485:2f} kJ/mol.")
    pass
