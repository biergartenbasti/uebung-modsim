# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os # Macht erstellen von Ordnern und Pfaden einfacher
import copy # Wird gebraucht um die Startstruktur unverändert zu lassen
from ase.visualize import view # import visualisation tools
from ase.io import read, write # Zum Lesen und Schreiben von Strukturdateien
from ase.constraints import ExpCellFilter # Wird für die Strukturoptimierung von periodischen Systemen benötigt
from ase.optimize import BFGS # Solver für die Berechnungen
from ase.md import Langevin  # NVT-Implementation für eine MD
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution # Wird fpr die Initialisierung der Geschwindigkeiten benötigt
from ase import units # Macht das Handling und die Umrechnung von Einheiten einfacher
from ase.io.trajectory import Trajectory # Wird benötigt um die Trajektorie zu speichern und an das Atoms object anzuhängen
from ase.io.vasp import write_vasp_xdatcar # Um Trajektorien zu XDATCAR zu schreiben
from mace.calculators import MACECalculator # Ein Rechner zur Berechnung
from mace.calculators import mace_mp # Das Foundation Model
import numpy as np
#import pandas as pd # vielleicht um Werte zu speicher
from pymatgen.io.vasp import Xdatcar


#%%
path_to_Li = 'Li/POSCAR_4x4x3_20vac'
path_to_Mg = 'Mg/POSCAR_4x4x3_20vac'
path_to_MgBulk = 'Mg/POSCAR_3x3x3Bulk'
path_to_LiPSCl = 'LiPSCl/POSCAR_2x2x2'
path_to_MgScSe =  'MgSc2Se4/POSCAR_defect'

structure_Li = read(path_to_Li,format='vasp')
structure_Mg = read(path_to_Mg,format='vasp')
structure_MgBulk = read(path_to_MgBulk,format='vasp')

structure_LiPSCl = read(path_to_LiPSCl,format='vasp')
structure_MgScSe = read(path_to_MgScSe,format='vasp')
#%% Calculator and relaxation
def relaxation_mace(structure, calc, fix_cell = True, f = 0.01):
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
    structure_data : DICTIONARY
        Beinhaltet die optimierte Struktur als ATOMS object und separat die berechnete Energie.
    '''
    # define calculator
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
    structure_data={'relaxed_structure':structure,'total_energy':structure.get_total_energy()}
    return structure_data

def MD_mace(structure, traj_out = 'output', T_start = 300.0, t_step = 3, n_steps = 1000, friction_coefficient = 1E-3):
    '''    
    Parameters
    ----------
    structure : ATOMS
        DESCRIPTION.
    traj_out : STRING, optional
        DESCRIPTION. The default is 'output'.
    T_start : FLOAT, optional
        DESCRIPTION. The default is 300.0.
    t_step : FLOAT, optional
        DESCRIPTION. The default is 3.
    n_steps : INT, optional
        DESCRIPTION. The default is 1000.
    friction_coefficient : FLOAT, optional
        DESCRIPTION. The default is 1E-3.

    Returns
    -------
    structure : ATOMS
        DESCRIPTION.
    '''
    
    # Geschwindigkeiten initialisieren
    MaxwellBoltzmannDistribution(structure, temperature_K=T_start)

    # Langevin-Setup für ein NVT-Ensemble
    dyn = Langevin(structure, timestep = t_step * units.fs, temperature_K = T_start, friction = friction_coefficient / units.fs)
    
    # Wir können einen Ordner anlegen, in dem die Trajektorien gespeichert werden
    try:
        os.makedirs(traj_out)
        print(f"Nested directories '{traj_out}' created successfully.")
    except FileExistsError:
        print(f"One or more directories in '{traj_out}' already exist.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{traj_out}'.")
    except Exception as e:
        print(f"An error occurred: {e}")
    # Definition einer Trajektorie, die angehängt werden kann
    traj = Trajectory(f'{traj_out}/{int(T_start)}K.traj', 'w', structure)
    dyn.attach(traj.write, interval=1)  # Writes to the trajectory file every 10 steps

    # Start der MD Simulation
    dyn.run(n_steps)
    traj.close()
    
    # Write als XDATCAR zusätzlich
    saved_traj = Trajectory(f'{traj_out}/{int(T_start)}K.traj')
    xdatcar_file = f'{traj_out}/XDATCAR_{int(T_start)}K'
    write_vasp_xdatcar(xdatcar_file, saved_traj, label=None)

    return structure
#%% Define calculators

mace_calc = mace_mp(model='medium-mpa-0',device='cuda',default_type='float32')
model_pt = "./MACE-matpes-r2scan-omat-ft.model"

mace_calc_bulk = MACECalculator(
    model_paths=model_pt,
    device='cuda',  # or 'cpu'
    default_dtype='float32'
)

#%% MD Simulation MgScSe Bulk
MgScSe_relaxed = relaxation_mace(structure_MgScSe, mace_calc, fix_cell = False, f = 0.01)

T_list = [800.0,1200.0,1600.0,2000.0]
traj_dir = 'MgSc2Se4Defect_Trajectories_3fs'
for Temp in T_list:
    print(f'Temperatur ist {Temp}')
    structure_tempo = copy.deepcopy(MgScSe_relaxed['relaxed_structure'])
    MD_mace(structure_tempo, traj_out = traj_dir, T_start = Temp, t_step = 3, n_steps = 8000, friction_coefficient = 1E-3)
    pass
#%% MD Simulation LiPSCl Bulk
LiPSCl_relaxed = relaxation_mace(structure_LiPSCl, mace_calc, fix_cell = False, f = 0.01)

T_list = [300.0,500.0,700.0,900.0,1100.0,1300.0]
traj_dir = 'LiPSCl_Trajectories'
for Temp in T_list:
    print(f'Temperatur ist {Temp}')
    structure_tempo = copy.deepcopy(LiPSCl_relaxed['relaxed_structure'])
    MD_mace(structure_tempo, traj_out = traj_dir, T_start = Temp, t_step = 2.5, n_steps = 1000, friction_coefficient = 1E-3)
    pass
#%% MD Simulation Lithium Surface
Li_relaxed=relaxation_mace(structure_Li, mace_calc, fix_cell = True, f = 0.01)

T_Li = [200.0,300.0,400.0,450.0]
Li_traj_dir = 'Li_Trajectories'
for Temp in T_Li:
    print(f'Temperatur ist {Temp}')
    structure_tempo = copy.deepcopy(Li_relaxed['relaxed_structure'])
    MD_mace(structure_tempo, traj_out = Li_traj_dir, T_start = Temp, t_step = 2.5, n_steps = 1000, friction_coefficient = 1E-3)
    pass

#%% MD Simulation Magnesium Surface
Mg_relaxed=relaxation_mace(structure_Mg, mace_calc, fix_cell = True, f = 0.01)

T_Mg = [200.0,300.0,400.0,500.0,600.0,700.0]
Mg_traj_dir = 'Mg_Trajectories'
for Temp in T_Mg:
    print(f'Temperatur ist {Temp}')
    structure_tempo = copy.deepcopy(Mg_relaxed['relaxed_structure'])
    MD_mace(structure_tempo, traj_out = Mg_traj_dir, T_start = Temp, t_step = 2.5, n_steps = 1000, friction_coefficient = 1E-3)
    pass
