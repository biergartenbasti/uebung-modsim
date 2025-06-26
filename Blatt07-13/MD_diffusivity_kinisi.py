#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:24:05 2024
Modified for ASE trajectory files

@author: sutz
"""

import numpy as np
import matplotlib.pyplot as plt
from kinisi.analyze import DiffusionAnalyzer
from kinisi.arrhenius import StandardArrhenius
import warnings
import os
import os.path
from ase.io import read
import corner
warnings.filterwarnings("ignore", category=UserWarning)

# import variables from bash parent
for key,value in os.environ.items():
    globals()[key]=value
    pass
MD_target='MgSc2Se4_Trajectories_1fs'
plot_title='MgSc2Se4 Bulk Diffusion 1 fs steps'

# define folder that contains trajectory files and create a list of them
trajectories_list = [f'{MD_target}/{el}' for el in os.listdir(MD_target) if el.endswith('.traj')]

np.random.seed(42)
rng = np.random.RandomState(42) # is not used right now because of reproducibility

p_params = {'specie': 'Mg',
          #'specie_indices': [97],
          'time_step': 1,
          'step_skip': 1,
          'progress': False,
          #'max_dt': 3000,
          #'min_dt': 100,
          #'n_steps': 200
          }
u_params = {'dimension':'xyz',
            'progress': False}
d_params={'progress': False}

time_skip=500.0

temps=np.array([]) # empty array, temps will be appended in the loop
D=[] # append diffusion coeffs
analyzers=[] #append analyzers

for traj_file in trajectories_list:
    traj_name = os.path.basename(traj_file) # get the base name of the trajectory file
    
    # Extract temperature from filename - you'll need to adjust this based on your naming convention
    # Assuming format like "simulation_temp_500K.traj" or similar
    # Please modify this line based on your actual file naming pattern:
    temp = float(traj_name.split('_')[-1].replace('K.traj', '').replace('.traj', ''))
    
    temps = np.append(temps, temp)
    
    # Read ASE trajectory
    trajectory = read(traj_file, ':')
    
    # Create DiffusionAnalyzer from ASE trajectory
    diff = DiffusionAnalyzer.from_ase(trajectory, parser_params=p_params, uncertainty_params=u_params)
    diff.diffusion(time_skip, diffusion_params={'random_state': None, 'progress': False})
    D.append(diff.D)
    analyzers.append(diff)

    # Create subplot figure for this temperature
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'{traj_name} - Temperature: {temp} K', fontsize=14)
    
    # Plot 1: MSD vs time (linear scale)
    axes[0, 0].errorbar(diff.dt, diff.msd, diff.msd_std)
    axes[0, 0].set_ylabel('MSD/Å$^2$')
    axes[0, 0].set_xlabel(r'$\Delta t$/ps')
    axes[0, 0].set_title('MSD vs Time (Linear)')
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: MSD vs time (log scale)
    axes[0, 1].errorbar(diff.dt, diff.msd, diff.msd_std)
    axes[0, 1].axvline(time_skip, color='g', label=f'time_skip={time_skip}')
    axes[0, 1].set_ylabel('MSD/Å$^2$')
    axes[0, 1].set_xlabel(r'$\Delta t$/ps')
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_xscale('log')
    axes[0, 1].set_title('MSD vs Time (Log-Log)')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()

    # Plot 3: MSD with credible intervals
    credible_intervals = [[15, 85], [2.5, 97.5], [0.15, 99.85]]
    alpha = [0.05, 0.1, 0.2]
    
    axes[1, 0].plot(diff.dt, diff.msd, 'k-', linewidth=2)
    for i, ci in enumerate(credible_intervals):
        axes[1, 0].fill_between(diff.dt,
                         *np.percentile(diff.distribution, ci, axis=1),
                         alpha=alpha[i],
                         color='#0173B2',
                         lw=0)
    axes[1, 0].set_ylabel('MSD/Å$^2$')
    axes[1, 0].set_xlabel(r'$\Delta t$/ps')
    axes[1, 0].set_title('MSD with Credible Intervals')
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Histogram of diffusion coefficients
    axes[1, 1].hist(diff.D.samples, density=True, bins=30, alpha=0.7)
    axes[1, 1].axvline(diff.D.n, c='k', linewidth=2, label=f'Mean D = {diff.D.n:.2e}')
    axes[1, 1].set_xlabel('$D$/cm$^2$s$^{-1}$')
    axes[1, 1].set_ylabel('$p(D$/cm$^2$s$^{-1})$')
    axes[1, 1].set_title('Diffusion Coefficient Distribution')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    print(f"{temp} K; D = {diff.D.n:.2e} ± {diff.D.ci()}")
    print(f"Intercept = {diff.intercept.n:.2e} ± {diff.intercept.ci()}")
    
    # Corner plot for this temperature
    corner.corner(diff._diff.flatchain, labels=['$D$/cm$^2$s$^{-1}$', 'intercept/Å$^2$'])
    plt.show()

#%%
# arrhenius parameters

s = StandardArrhenius(temps, D, bounds=((1E-3,1),(1e-20,1E-1)))
s.max_likelihood('mini')
s.mcmc()

#%% visualise probability of parameter distribution
plt.close('all')
unit_cms=' cm\u00b2 s\u207B\u00b9'
fig=corner.corner(s.flatchain, labels=['$E_a$/eV', '$A$/cm$^2$s$^{-1}$'])
fig.set_size_inches(8,8) # dont know how to change the size otherwise
corner.overplot_lines(fig, s.variable_medians, color="orange")
corner.overplot_points(fig, s.variable_medians[None], marker="s", color="orange")
corner.overplot_lines(fig, s.variable_modes, color="green")
corner.overplot_points(fig, s.variable_modes[None], marker="s", color="green")

# place a text box in upper left in axes coords
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
inplot_text="Median Ea = "+np.format_float_scientific(np.round(s.variable_medians[0],3))+" eV\nMedian A = "+np.format_float_scientific(s.variable_medians[1],precision=2)+unit_cms
inplot_text2="\nMode Ea = "+np.format_float_scientific(np.round(s.variable_modes[0],3))+" eV\nMode A = "+np.format_float_scientific(s.variable_modes[1],precision=2)+unit_cms
fig.text(0.6,0.75, inplot_text+inplot_text2,fontsize=14,verticalalignment='center',ha='left',bbox=props)
fig.savefig(MD_target+'/ArrheniusParamsDistribution_kinisi'+plot_title+'.png',dpi=1000)

# print arrhenius parameters information and save to file
print("Modes are green, medians are orange")
print("Modes are Ea / eV: "+np.format_float_scientific(s.variable_modes[0])+", A / cm2 s-1: "+np.format_float_scientific(s.variable_modes[1],precision=2))
print("Medians are Ea / eV: "+np.format_float_scientific(s.variable_medians[0])+", A / cm2 s-1: "+np.format_float_scientific(s.variable_medians[1],precision=2))
# diffusion coefficient at 300K
D300K=np.format_float_scientific(s.variable_medians[1]*np.exp(-s.variable_medians[0]/300/1E-5/8.62),precision=2)
print("Using median values, D / cm2 s-1 at 300K is: "+D300K)
print("Median Absolute Deviation of A / cm2 s-1 is: "+np.format_float_scientific(np.median(np.absolute(s.variables[0] - np.median(s.variables[0]))),precision=2))
print("Median Absolute Deviation of Ea / eV is: "+np.format_float_scientific(np.median(np.absolute(s.variables[1] - np.median(s.variables[1]))),precision=2))
with open(MD_target+'/MD_Parameters_kinisi.txt', 'w') as f:
    f.write("Modes are green, medians are orange.\n")
    f.write("Modes are Ea / eV: "+np.format_float_scientific(s.variable_modes[0])+", A / cm2 s-1: "+np.format_float_scientific(s.variable_modes[1])+'\n')
    f.write("Medians are Ea / eV: "+np.format_float_scientific(s.variable_medians[0])+", A / cm2 s-1: "+np.format_float_scientific(s.variable_medians[1])+'\n')
    f.write("Median Absolute Deviation of A / cm2 s-1 is: "+np.format_float_scientific(np.median(np.absolute(s.variables[0] - np.median(s.variables[0]))))+'\n')
    f.write("Median Absolute Deviation of Ea / eV is: "+np.format_float_scientific(np.median(np.absolute(s.variables[1] - np.median(s.variables[1]))))+'\n')

    f.write("\nUsing median values, D / cm2 s-1 at 300K is: "+D300K)
    pass
f.close()

# plot arrhenius function
fig2,ax2=plt.subplots(figsize=(10,8))

# plot data points
plt.errorbar(1000/s.x, s.y.n, s.y.ci(), marker='o', ls='', color='k', zorder=10)

Tx = np.linspace(300, 2500, 200) # temps
Dy = s.variable_medians[1]*np.exp(-s.variable_medians[0]/Tx/1E-5/8.62) # Ds
Dy2=s.variable_medians[0]/8.62/1E5*1000/Tx
plt.plot(1000/Tx, Dy,color='red')
plt.yscale('log')
plt.xlabel('$1000T^{-1}$/K$^{-1}$')
plt.ylabel('$D$/cm$^2$s$^{-1}$')
plt.title(plot_title)

# place a text box in upper left in axes coords
ax2.text(0.75,0.95, inplot_text+'\nD$_{300K}$ = '+D300K+unit_cms, transform=ax2.transAxes, fontsize=18,
        verticalalignment='top',ha='center',bbox=props)
plt.savefig(MD_target+'/ArrheniusPlot_kinisi'+plot_title+'.png',dpi=1000)
plt.show()
#%%
# Add this before the Arrhenius fitting to examine your data
print("Temperature and D values:")
for i, (temp, d_val) in enumerate(zip(temps, D)):
    ci_lower, ci_upper = d_val.ci()
    uncertainty = (ci_upper - ci_lower) / 2  # Half-width of confidence interval
    print(f"T = {temp} K: D = {d_val.n:.2e} ± {uncertainty:.2e} cm²/s")
    print(f"    CI: [{ci_lower:.2e}, {ci_upper:.2e}]")

# Check for outliers or unrealistic values
plt.figure(figsize=(10, 6))
d_values = [d.n for d in D]
d_errors = [(d.ci()[1] - d.ci()[0])/2 for d in D]  # Half-width of CI as error

plt.errorbar(temps, d_values, d_errors, 
             marker='o', capsize=5, capthick=2)
plt.xlabel('Temperature (K)')
plt.ylabel('D (cm²/s)')
plt.yscale('log')
plt.title('Diffusion Coefficients vs Temperature')
plt.grid(True, alpha=0.3)
plt.show()