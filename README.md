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

### Simulation Workflow

The simulation proceeds as follows:

**1. Lattice Generation**
   * Construction of an FCC crystal
   * Periodic boundary conditions (PBC) applied
   * Computation of (initial) pairwise scalar distances and distance vectors
2. Force and Energy Calculation
   * Lennard–Jones forces computed for all atom pairs within the cutoff radius
   * Potential energy calculated per atom and for the full lattice

Under periodic boundary conditions, all atoms exhibit similar potential energy due to the absence of surface effects.

3. Thermodynamic Initialization
   * Random velocity distribution generated
   * Center-of-mass velocity removed
   * Velocities rescaled to match a target temperature $$T_0 = 300 K$$

This is an initialization step only, no thermostat is applied during dynamics

4. Time Integration (Velocity-Verlet)

The system evolves according to:

$$
\mathbf{r}_{t+\Delta t} =
\mathbf{r}_t + \Delta t \, \mathbf{v}_{t+\frac{1}{2}\Delta t}
$$

$$
\mathbf{v}_{t+\frac{1}{2}\Delta t} =
\mathbf{v}_t + \frac{1}{2}\Delta t \, \mathbf{a}_t
$$

with

$$
\mathbf{a}_t = \frac{\mathbf{F}_t}{m}
$$

Positions, velocities, energies, and temperature are recorded at each timestep.



