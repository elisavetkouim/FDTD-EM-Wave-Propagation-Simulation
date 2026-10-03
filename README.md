# FDTD-EM-Wave-Propagation-Simulation

*This project was part of the course Computational Electromagnetics at Aristotle University of Thessaloniki.*
---
## Overview
The project inculdes MATLAB implementations of the Finite Differences in the Time Domain (FDTD) method for Electromagnetic wave propagation for two different cases with PML boundary conditions.

More specifically:

### 1. **2D TM ($E_z, H_x, H_y$) Wave Propagation in Parallel-Plate Waveguide**
  - 2D Yee grid formulation for the three field compenents
  - Time updating loop for the field components based on the Discretized Maxwell Equations, derived from implementation of central differences in time and space, on Maxwell's equations
  - Hard source implemenation in the center of the computational space
  - Evaluation of the Courant stability condition demonstrating field divergence for $\Delta t > \Delta t_c$
  - Implementation of Perfectly Matched Layer (PML) boundary conditions to absorb outgoing waves and eliminate boundary reflections
  - Computation of phase velocity ($v_p$) for different grid resolutions ($\lambda/5, \lambda/10, \lambda/20$) and comparison against analytical value
<img width="400" height="137" alt="FDTD_waveguide_animation3" src="https://github.com/user-attachments/assets/8055e1c1-f94b-4f05-9d3b-7c1d99809651" />

### 2. **1D TM ($E_z, H_x$) Wave Propagation & interaction with time-varying dielectric medium**
  - Computation of the Modified Discretized Maxwell Equations for time-dependent perimittivity, derived from implementation of central differences in time and space, on Maxwell's equations
  - 1D Yee grid formulation for the two field compenents
  - Time updating loop for the field components
  - Hard source implemenation on the left of the dielectric plate with time-varying permittivity
  - PML implementation on the left and right boundaries of the computational space, study of three different cases of layer thickness (8, 16 & 100 cells)
<img width="400" height="137" alt="fdtd_simulation_PML_16" src="https://github.com/user-attachments/assets/faa42ab2-1348-443c-8239-935b9ef669d7" />

---
*The project report is in Greek*
