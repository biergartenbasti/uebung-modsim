# ��� MD Simulation Interactive Lab

[![Launch Simulation](https://img.shields.io/badge/Launch-Interactive%20Simulation-blue?style=for-the-badge&logo=jupyter)](https://mybinder.org/v2/gh/biergartenbasti/uebung-modsim/HEAD?urlpath=voila%2Frender%2Fnotebooks%2FMD_Simulation.ipynb)
[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-green?style=for-the-badge&logo=github)](https://biergartenbasti.github.io/uebung-modsim/)

## ��� Quick Start

**[��� Launch Interactive Simulation](https://biergartenbasti.github.io/uebung-modsim/)** - Beautiful landing page with direct access

**[��� Direct Binder Link](https://mybinder.org/v2/gh/biergartenbasti/uebung-modsim/HEAD?urlpath=voila%2Frender%2Fnotebooks%2FMD_simulation.ipynb)** - Skip straight to the simulation

## ��� About

This interactive molecular dynamics simulation allows you to:

- **Explore different integration methods**: Verlet, Velocity Verlet, Langevin, Monte Carlo
- **Adjust simulation parameters in real-time** with intuitive sliders
- **Visualize particle trajectories** and energy evolution
- **Analyze diffusion properties** with distance-squared plots
- **No installation required** - runs entirely in your browser!

## ��� Interactive Controls

| Parameter | Range | Description |
|-----------|--------|-------------|
| **Integrator** | 0-3 | Integration method (0=Verlet, 1=Velocity Verlet, 2=Langevin, 3=Monte Carlo) |
| **Temperature** | 200-500 K | Simulation temperature |
| **Gamma** | 10-15 | Friction coefficient for Langevin dynamics |
| **Steps** | 50-500 | Number of simulation timesteps |
| **E_kin** | 10-50 meV | Initial kinetic energy |
| **Angle** | 0-90° | Initial velocity direction |

## ���️ Technical Details

- **Built with**: Jupyter, Voilà, ipywidgets, matplotlib, numpy, scipy
- **Hosted on**: GitHub Pages + Binder
- **No server required**: Runs client-side after initial setup
- **Open source**: Full code available in this repository

## ��� Repository Structure
