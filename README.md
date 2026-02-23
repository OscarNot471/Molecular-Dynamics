# Molecular-Dynamics
## Lennard–Jones FCC Crystal with Velocity-Verlet Integration

<p align="center">
  <img src="results/Molecular_dynamics_animation_short.gif" width="750">
</p>

## Project Overview
This project was developed as part of the course “Introduction to Modelling in Physics”, which consists of 14 computational modelling assignments focused on simulating physical systems and solving differential equations.

In practical sessions 6 and 7, the objective was to implement from scratch a molecular dynamics (MD) simulation of a crystalline solid using:
* Pairwise Lennard–Jones potential
* Thermodynamic initialization of velocities
* Time integration via the Velocity-Verlet algorithm
* Energy and temperature monitoring
* Trajectory export for visualization in OVITO
The implementation is fully vectorized (except for the time loop) and does not rely on external MD libraries.
