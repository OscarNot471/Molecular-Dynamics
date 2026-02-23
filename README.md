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

## Physical Model

We simulate a 3D face-centered cubic (FCC) lattice of copper atoms.

### Interaction Potential

Atoms interact through the classical Lennard–Jones pair potential:

$$ V(r) = 4\epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^6 \right] $$

Forces are computed as the analytical gradient of the potential.

A cutoff radius is applied to reduce computational cost.







