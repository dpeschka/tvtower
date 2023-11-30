# Elastic Fernsehturm.
Simple demonstration of damped Hamiltonian elastodynamics with FEniCS.

[Animated TV tower from elastic Fernsehturm code](tvtower.gif)

Author: Dirk Peschka

## What does this code do?

1. Creation of a simplified two-dimensional representation of the Berlin Fernsehturm tower as a finite element mesh.
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

where the Frechet derivative in FEniCS is conveniently obtained via `derivative(H,q,dq)`. In the actual implementation we use the displacement $w(x)$ instead of the deformation $\chi(x):=x+w(x)$. Then the elastodynamic problem is solved via

$$
a(\eta,\xi) - b(\xi,\partial_t q)=0,
$$

where $a(\eta,xi)=j(\eta,\xi)-k(\eta,\xi)$ consists of a skew-symmetric contribution $j$ that covers reversible (Hamiltonian) dynamics and a symmetric positive contribution $k$ that covers irreversible (Onsager) dynamics. We use 

$$
j(\eta,\xi)=\int_\Omega \eta_p\cdot\xi_\chi - \xi_p\cdot\eta_\chi\,\mathrm{d}x,
$$

and

$$
k(\eta,\xi)=\int_\Omega \mu\nabla\eta_p\cdot\nabla \xi_p \,\mathrm{d}x.
$$

We discretize via P2 FEM for all function and a Crank-Nicolson scheme in time.

