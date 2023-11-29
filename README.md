# TV tower of Berlin.
Simple demonstration of damped Hamiltonian elastodynamics with FEniCS.

## What does this code do?

1. Creation of a simplified two-dimensional representation of the Berlin TV tower as a finite element mesh.
2. Somewhat compact representation of a structure-preserving damped-Hamiltonian system of elastodynamics.
3. Solve and visualize the solution using the finite element library FEniCS (legacy version 2019.1.0)

## What do you need to run this code?

- Python
- FEniCS legacy 2019.1.0 with mshr
- Matplotlib

## Short explanation of discretization

The starting point of this discretization is a formulation of elastodynamics using momentum $p$ and deformation $\chi$.
