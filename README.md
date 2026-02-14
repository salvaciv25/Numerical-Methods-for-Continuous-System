# Finite Element Solvers for Advection-Diffusion and Stokes Problems

This repository contains MATLAB implementations for the numerical solution of partial differential equations (PDEs) related to fluid dynamics. The project is divided into two main parts: a **2D Convection-Diffusion solver** with upwind stabilization and a **Stokes flow solver** using different stable/stabilized FEM elements.

## Features

### 1. Convection-Diffusion Solver (Project 1)

Solves the steady-state Advection-Diffusion equation on a unit square domain .

* **Method:** Galerkin Finite Element Method (FEM) with  linear triangular elements.
* **Stabilization:** Streamline Upwind Diffusion (SUD) to handle high Peclet numbers and prevent numerical oscillations.
* **Boundary Conditions:** Dirichlet conditions handled via boundary lifting operators.
* **Scenarios:** Tested across different regimes, from pure diffusion to convection-dominated flows ().

### 2. Stokes Equation Solver (Project 2)

Solves the incompressible Stokes equations for velocity  and pressure .

* **Algorithms Implemented:**
1. **P1/P1:** Standard linear elements for both velocity and pressure (unstable, used for comparison).
2. **P1/P1 GLS (Galerkin Least Squares):** Stabilized formulation using a pressure-stabilizing term ().
3. **Mini-Element (Bubble Scheme):** Enrichment of the velocity space with cubic bubble functions to satisfy the Info-Sup (LBB) condition.


* **Test Case:** Validated on the **Lid-Driven Cavity** problem.
* **Error Analysis:** Comparison of residual errors across multiple mesh refinements.

---

## Implementation Details

The code is modularized into specialized functions to ensure clarity and reusability:

| Function | Description |
| --- | --- |
| `input_var.m` | Parses mesh data and connectivity from the `/mesh` folder. |
| `LocalB.m` | Evaluates  basis functions, gradients, and element areas. |
| `MatrixB.m` | Builds global Stiffness (), Transport (), and Stabilization () matrices. |
| `lifting.m` | Applies the boundary lifting operator for non-homogeneous Dirichlet BCs. |
| `stokes_solver.m` | Unified function for P1/P1, GLS, and Bubble schemes. |

---

## Key Results

### Convection-Diffusion Stability

When the diffusion coefficient  is small, the standard Galerkin method produces spurious oscillations. The **SUD stabilization** effectively reduces these artifacts, allowing for accurate solutions even at high Peclet numbers.

### Stokes Convergence

The **P1/P1 GLS** and **Bubble** schemes provide stable pressure fields, unlike the standard P1/P1 approach which suffers from pressure instability.

* **GLS Tuning:** Observations show that  is optimal for the Lid-Driven Cavity; higher values lead to excessive smoothing of the solution.

---

## How to Run

1. Clone the repository.
2. Add the mesh folders to your MATLAB path.
3. Execute the main scripts:
* `fem_solver.m` for Convection-Diffusion analysis.
* `Project2_CS.m` for Stokes flow simulations.
* `Project2_CS_LidCavityPb.m` for Lid Cavity problem.


4. The results will automatically display plots for velocity components (), pressure (), and error convergence.
