# 2D Ising Model Monte Carlo Simulation

A Fortran 95 implementation of a Monte Carlo simulation of the two-dimensional Ising model using the Metropolis algorithm.

## Overview

The simulation models a square lattice of spins with periodic boundary conditions and studies its thermodynamic behaviour as a function of temperature.

The simulation calculates:
- Magnetization and absolute magnetization
- Energy
- Magnetic susceptibility
- Specific heat
- Binder cumulant

## Method

- 2D square lattice (16 × 16)
- Metropolis Monte Carlo algorithm
- Periodic boundary conditions
- Thermalization before measurements
- Temperature-dependent sampling of thermodynamic observables

## Implementation

The code is written in **Fortran 95** and uses modular programming with separate modules for simulation constants and lattice operations.

## Output

The simulation writes the temperature-dependent observables to a `.dat` file for further analysis and plotting.
