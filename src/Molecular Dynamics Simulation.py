import numpy as np
import matplotlib.pyplot as plt
import sys


def generate_fcc_lattice(nx, ny, nz, a):
    """
    Generate a three-dimensional Face-Centered Cubic (FCC) lattice.

    The lattice is constructed by replicating the FCC unit cell along the
    x, y, and z directions. For each unit cell of lattice parameter 'a',
    four atoms are added: one at the corner and three at the face-centered
    positions.

    Returns an array of shape (N, 3) with the Cartesian coordinates of all atoms.
    """
    positions = []

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # Corner atom
                positions.append([i * a, j * a, k * a])

                # Face-centered atoms
                positions.append([i * a + 0.5 * a, j * a + 0.5 * a, k * a])
                positions.append([i * a + 0.5 * a, j * a, k * a + 0.5 * a])
                positions.append([i * a, j * a + 0.5 * a, k * a + 0.5 * a])

    return np.array(positions)

"""
Compute the Lennard-Jones potential for a given interparticle distance r.
"""

def lennard_jones_potential(r, sigma=2.3151, epsilon=0.167, n=12, m=6):

    # Avoid division by zero
    r_inv = np.copy(r)
    r_inv[r != 0] = sigma / r_inv[r != 0]

    return 4 * epsilon * (r_inv**n - r_inv**m)

"""
Compute the interparticle forces from the derivative of the Lennard-Jones potential.
"""

def lennard_jones_force(distances, displacement_vectors, cutoff_radius,
                        sigma=2.3151, epsilon=0.167, n=12, m=6):

    # Avoid division by zero
    r_inv = np.zeros_like(distances)
    r_inv[distances != 0] = 1 / distances[distances != 0]

    # Compute unit vectors between particles
    unit_vectors = np.zeros_like(displacement_vectors)
    norm = np.ones_like(distances)
    norm[distances != 0] = distances[distances != 0]
    unit_vectors = displacement_vectors / norm[..., np.newaxis]

    # Force magnitude from LJ derivative
    force_magnitude = -4 * epsilon * (
        n * (sigma**n) * (r_inv**(n+1)) -
        m * (sigma**m) * (r_inv**(m+1))
    )

    # Pairwise force vectors
    pair_forces = force_magnitude[..., np.newaxis] * unit_vectors

    # Apply cutoff radius and remove self-interaction
    pair_forces[(distances > cutoff_radius) | (distances == 0)] = 0.0

    # Total force on each particle
    total_force = np.sum(pair_forces, axis=1)

    return total_force

"""
Compute pairwise distances and displacement vectors for a given lattice
(with free boundary conditions).
"""

def compute_free_boundary_distances(positions):

    num_particles = len(positions)
    distances = np.zeros((num_particles, num_particles))
    displacement_vectors = np.zeros((num_particles, num_particles, 3))

    for i in range(num_particles):
        r = positions - positions[i]

        distances[i] = np.linalg.norm(r, axis=1)
        displacement_vectors[i, :, :] = r[:, :]

    return distances, displacement_vectors


"""
Compute pairwise distances and displacement vectors using
periodic boundary conditions (minimum image convention).
"""

def compute_periodic_distances(positions, nx, ny, nz, a):

    num_particles = len(positions)

    Lx = nx * a
    Ly = ny * a
    Lz = nz * a

    distances = np.zeros((num_particles, num_particles))
    displacement_vectors = np.zeros((num_particles, num_particles, 3))

    for i in range(num_particles):

        r_periodic = positions - positions[i]

        # Apply minimum image convention
        r_periodic[:, 0] -= Lx * np.round(r_periodic[:, 0] / Lx)
        r_periodic[:, 1] -= Ly * np.round(r_periodic[:, 1] / Ly)
        r_periodic[:, 2] -= Lz * np.round(r_periodic[:, 2] / Lz)

        displacement_vectors[i, :, :] = r_periodic[:, :]
        distances[i] = np.linalg.norm(r_periodic, axis=1)

    return distances, displacement_vectors


"""
Given the pairwise distances and displacement vectors of a configuration,
compute the Lennard-Jones potential energy and the forces acting on each atom,
taking into account a cutoff radius.
"""

def compute_energy_and_forces(distances, displacement_vectors,
                              sigma=2.3151, epsilon=0.167, n=12, m=6):

    cutoff_radius = 4 * sigma

    mask = (distances > 0) & (distances < cutoff_radius)

    # Pairwise potential energy
    energy_matrix = lennard_jones_potential(distances, sigma, epsilon, n, m)

    # Remove self-interaction and apply cutoff
    energy_matrix = np.where(mask, energy_matrix, 0)

    # Total potential energy
    total_energy = np.sum(energy_matrix) / 2.0

    # Energy per particle
    energy_per_particle = np.sum(energy_matrix, axis=1) / 2.0
    num_particles = distances.shape[0]
    mean_energy = total_energy / num_particles

    # Forces (already include cutoff and r != 0 handling)
    forces = lennard_jones_force(distances, displacement_vectors,
                                 cutoff_radius, sigma, epsilon, n, m)

    return total_energy, mean_energy, energy_per_particle, forces


"""
Plot the 3D atomic configuration using a colormap to represent the per-particle energy.
"""

def plot_3d_energy_colormap(positions, energy_per_particle, point_size=3):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    p = ax.scatter(
        positions[:, 0], positions[:, 1], positions[:, 2],
        s=120, c=energy_per_particle, cmap=plt.cm.viridis
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    fig.colorbar(p, ax=ax)
    plt.show()


"""
Plot the force vectors acting on each atom as a 3D quiver plot.
"""

def plot_3d_forces(positions, forces, colorbar_label=r"$F$"):

    vectors = -forces

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("Force on each atom")

    p = ax.quiver(
        positions[:, 0], positions[:, 1], positions[:, 2],
        vectors[:, 0], vectors[:, 1], vectors[:, 2],
        length=1, normalize=True
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    fig.colorbar(p, ax=ax, label=colorbar_label)
    plt.show()

    return


"""
Thermodynamic initialization
"""

# Kinetic energy
def compute_kinetic_energy(velocities, mass):
    return np.sum(0.5 * mass * (np.linalg.norm(velocities, axis=1)**2))


# Temperature
def compute_temperature(velocities, mass, kb):
    return (2/3) * compute_kinetic_energy(velocities, mass) / (kb * (4 * len(velocities)))


"""
Generate a random velocity distribution with a target temperature.
The center-of-mass velocity is removed and the distribution is rescaled
to match the desired temperature.
"""

def initialize_velocities(num_particles, target_temperature, mass, kb):

    velocities = np.random.random_sample((num_particles, 3)) - 0.5

    # Remove center-of-mass velocity
    v_cm = np.sum(velocities, axis=0) / num_particles
    velocities -= v_cm

    # Temperature of the unscaled distribution
    random_temperature = compute_temperature(velocities, mass, kb)

    # Rescale velocities to match target temperature
    scale = np.sqrt(target_temperature / random_temperature)
    velocities *= scale

    # Optional check
    test_temperature = compute_temperature(velocities, mass, kb)
    print(test_temperature)

    return velocities


"""
Verlet integration algorithm.

Time evolution of the atomic system. Positions, velocities,
forces, potential energy, kinetic energy and temperature
are computed at each time step.
"""

def verlet(positions, initial_velocities, num_steps, time_step, mass):

    num_particles = len(positions)

    # Initial distances and displacement vectors (periodic boundaries)
    distances, displacement_vectors = compute_periodic_distances(
        positions, nx, ny, nz, a
    )

    # Position history
    r = np.zeros((num_steps, num_particles, 3))
    r[0] = positions

    # Velocity history
    v = np.zeros((num_steps, num_particles, 3))
    v[0] = initial_velocities

    # Initial forces and energies
    total_energy, mean_energy, energy_per_particle, forces = \
        compute_energy_and_forces(distances, displacement_vectors)

    # First half-step velocities
    v_half = initial_velocities + 0.5 * time_step * (forces / mass)

    # Energy and temperature arrays
    potential_energy = np.zeros(num_steps)
    kinetic_energy = np.zeros(num_steps)
    temperature = np.zeros(num_steps)

    potential_energy[0] = total_energy
    kinetic_energy[0] = compute_kinetic_energy(v[0], mass)
    temperature[0] = compute_temperature(v[0], mass, kb)

    # Time evolution loop
    for i in range(1, num_steps):

        # Progress indicator (update every 1%)
        if i % (num_steps // 100 or 1) == 0:
            progress = int(100 * i / num_steps)
            sys.stdout.write(f"\rMD simulation: {progress}% completed")
            sys.stdout.flush()

        # Update positions
        r[i] = r[i-1] + time_step * v_half

        # Recompute distances and forces
        distances, displacement_vectors = compute_periodic_distances(
            r[i], nx, ny, nz, a
        )

        total_energy, mean_energy, energy_per_particle, forces = \
            compute_energy_and_forces(distances, displacement_vectors)

        # Update velocities
        v[i] = v_half + 0.5 * time_step * (forces / mass)

        # Energies and temperature
        potential_energy[i] = total_energy
        kinetic_energy[i] = compute_kinetic_energy(v[i], mass)
        temperature[i] = compute_temperature(v[i], mass, kb)

        # Next half-step
        v_half = v[i] + 0.5 * time_step * (forces / mass)

    sys.stdout.write("\rMD simulation: 100% completed\n")

    return r, v, potential_energy, kinetic_energy, temperature


############################################
# MAIN
############################################

# Lattice parameters
a = 3.603
N = 4
nx = N
ny = N
nz = N

atomic_mass_uma = 63.55              # Copper atomic mass (amu)
uma_to_evfsa2 = 103.64               # Conversion factor: amu → eV·fs^2 / Å^2
mass = atomic_mass_uma * uma_to_evfsa2

kb = 8.617e-5                        # Boltzmann constant (eV/K)

initial_temperature = 300            # Initial temperature (K)

num_steps = 1000                     # Total number of time steps
time_step = 1                        # Time step (fs)

# Lattice construction
positions = generate_fcc_lattice(nx, ny, nz, a)

# Initial distances and energies
distances, displacement_vectors = compute_periodic_distances(
    positions, nx, ny, nz, a
)

total_energy, mean_energy, energy_per_particle, forces = \
    compute_energy_and_forces(distances, displacement_vectors)

print("Total lattice energy:", total_energy, "eV")
print("Mean energy per atom:", mean_energy, "eV/atom")

plot_3d_energy_colormap(positions, energy_per_particle)
plot_3d_forces(positions, forces)

# Initialize velocities at target temperature
initial_velocities = initialize_velocities(
    len(positions), initial_temperature, mass, kb
)

# Time evolution using Verlet integration
r, v, potential_energy, kinetic_energy, temperature = \
    verlet(positions, initial_velocities, num_steps, time_step, mass)

# Shift potential energy reference
potential_energy_shifted = potential_energy - potential_energy[0]

# Total energy over time
total_energy_time = kinetic_energy + potential_energy_shifted

# Time array
time = np.linspace(0, num_steps * time_step, num_steps)


plt.figure()
plt.title("Lattice energies as a function of time")
plt.plot(time, total_energy_time, color="r", label="Total energy")
plt.plot(time, kinetic_energy, color="b", label="Kinetic energy")
plt.plot(time, potential_energy_shifted, color="g", label="Potential energy")
plt.xlabel("Time (fs)")
plt.ylabel("Energy (eV)")
plt.legend(loc="best")
plt.show()

plt.figure()
plt.title("Temperature as a function of time")
plt.plot(time, temperature)
plt.xlabel("Time (fs)")
plt.ylabel("Temperature (K)")
plt.show()


def save_xyz_trajectory(positions_time, time, num_particles, filename="trajectory"):

    positions_time = np.round(positions_time, 5)

    with open(filename + ".xyz", "w") as output:

        for i in range(len(time) // 2):

            output.write(f"{num_particles}\n")
            output.write(f"Atoms. Timestep: {int(time[2*i])}\n")

            for j in range(num_particles):
                output.write(
                    f"Cu  {positions_time[i][j][0]}  "
                    f"{positions_time[i][j][1]}  "
                    f"{positions_time[i][j][2]}\n"
                )

    return


save_xyz_trajectory(r, time, len(positions))