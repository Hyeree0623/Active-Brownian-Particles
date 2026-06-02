# Active Brownian Particle (ABP) Simulation

A high-performance C implementation of a two-dimensional Active Brownian Particle (ABP) simulation using overdamped Brownian dynamics.

The code simulates self-propelled particles interacting through a purely repulsive Weeks-Chandler-Andersen (WCA) potential under periodic boundary conditions.

---

## Features

- 2D Active Brownian Particle dynamics
- Rotational diffusion of particle orientation
- Self-propulsion with constant speed
- Weeks-Chandler-Andersen (WCA) interactions
- Periodic boundary conditions (PBC)
- Cell-list neighbor search for efficient force calculations
- Unwrapped and wrapped particle coordinates
- GRO trajectory output for visualization and analysis
- Scalable to several thousand particles

---

## Physical Model

Each particle evolves according to

dx/dt = μF + v₀n + η

where

- μF : interaction force term
- v₀ : self-propulsion velocity
- n : particle orientation vector
- η : translational Brownian noise

The particle orientation evolves through rotational diffusion

dθ/dt = √(2Dᵣ) ξ(t)

where

- Dᵣ : rotational diffusion coefficient
- ξ(t) : Gaussian white noise

---

## Interaction Potential

Particles interact through the Weeks-Chandler-Andersen (WCA) potential

U(r) =
4ε[(σ/r)¹² − (σ/r)⁶] + ε

for

r < 2^(1/6)σ

and

U(r) = 0

otherwise.

This provides purely repulsive excluded-volume interactions.

---

## Simulation Parameters

Default parameters:

| Parameter | Value |
|------------|--------|
| Number of particles | 5000 |
| Box size | 100 × 100 |
| Time step | 1×10⁻⁵ |
| Total steps | 10,000,000 |
| Self-propulsion speed (v₀) | 5.0 |
| Translational diffusion | 0.05 |
| Rotational diffusion | 0.15 |
| WCA ε | 1.0 |
| WCA σ | 1.0 |

---

## Algorithm Workflow

Initialization
    ↓
Random particle placement
    ↓
Build cell list
    ↓
Compute WCA forces
    ↓
Update particle orientation
    ↓
Apply active propulsion
    ↓
Apply Brownian motion
    ↓
Periodic boundary conditions
    ↓
Update cell list
    ↓
Save trajectory
    ↓
Repeat

---

## Output

### Trajectory File

The simulation produces:

coord.gro

containing

- particle coordinates
- simulation time
- box dimensions

which can be visualized using:

- VMD
- OVITO
- MDAnalysis
- custom Python scripts

---

## Performance

The code uses a cell-list neighbor searching algorithm.

Without cell lists:

O(N²)

With cell lists:

O(N)

for sufficiently large systems.

This allows efficient simulations of several thousand particles.

---

## Compilation

Using GCC:

```bash
gcc -O3 -lm ABP.c -o ABP
