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

The starting point of this discretization is a formulation of elastodynamics using momentum $p:\Omega\to\mathbb{R}^d$ and deformation $\chi:\Omega\to\mathbb{R}^d$, where $\Omega\subset\mathbb{R}^d$ has the shape of the TV tower and $d=2$. We put these functions into a composite vector $q=(\chi,p)$ that depends on time, i.e. $q=q(t)$ for $0\le t\le T$ with initial data $q_0=q(t=0)$. The discretization of the damped Hamiltonian system starts from a Hamiltonian (free energy)

$$
\mathcal{H}(q)=\int_\Omega \frac{|p|^2}{2\varrho} + W(F)\,\mathrm{d}x
$$

using the deformation gradient $F=\nabla\chi$. For this we define extra variables $\mathbf{\eta}=(\eta_\chi,\eta_p)$ and identify forces via 

$$
b(\eta,v)=\langle \mathrm{D}\mathcal{H}(q),v\rangle
$$

where the Frechet derivative in FEniCS is conveniently obtained via `derivative(H,q,dq)`.
