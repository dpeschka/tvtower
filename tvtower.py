from fenics import *
from pylab import plt
import numpy as np
set_log_level(50)

output_folder = "./output/"
file1 = File(output_folder + "displacement.pvd")
file2 = File(output_folder + "momentum.pvd")
file3 = File(output_folder + "velocity.pvd")
file4 = File(output_folder + "force.pvd")

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True,"eliminate_zeros": True,
               "precompute_basis_const": True,"precompute_ip_const": True}  

# Output routine
def output(q,t=0):
    # Output
    fw,fp,w,p = q.split()
    fw.rename('fw','Force')
    fp.rename('fp','Velocity')
    w.rename('w', 'Displacement')
    p.rename('p', 'Momentum')

    tmp  = Function(W)
    ALE.move(mesh, q.sub(2))
    file1 << (w, t)
    file2 << (p, t)
    file3 << (fp, t)
    file4 << (fw, t)
    tmp.assign(-q)
    ALE.move(mesh, tmp.sub(2))

# Material constants
mu        = Constant(500.0) # neo-Hookean elastic modulus
viscosity = Constant(0.05)  # viscosity
T         = 10.0            # final time
n_steps   = 100             # time steps
dt        = T/(n_steps)     # time step size

# Define mesh & function spaces
mesh = Mesh("tvtower.xml")
P2   = VectorElement("P", mesh.ufl_cell(), 2)
W    = FunctionSpace(mesh, MixedElement([P2,P2,P2,P2]))

# Boundary conditions, forcing terms, and initial conditions
initial = Expression(("0","0","0","0","0", "0", "x[1]*c", "0"),c=0.15,degree=2)
bc1 = DirichletBC(W.sub(2),Constant((0,0)),'on_boundary && near(x[1], 0, 1e-8)')
bc2 = DirichletBC(W.sub(3),Constant((0,0)),'on_boundary && near(x[1], 0, 1e-8)')

# Solve single time step
def evolve(old_q, dt):
    # unknowns, test functions and previous time step
    q,dq = Function(W),TestFunction(W)
    fw,fp,w,p                  = split(q)
    dfw,dfp,dw,dp              = split(dq)
    old_fw,old_fp,old_w, old_p = split(old_q)
    
    # Hamiltonian = kinetic + potential energy
    F = Identity(2) + grad(w)
    e_kinetic     = 0.5*p**2
    e_potential   = (mu/2)*(tr(F.T*F-Identity(2)) - 2*ln(det(F)))
    H             = (e_kinetic + e_potential)*dx
    
    # damped Hamiltonian formulation with dual variables fw,fp
    Res  = inner( (w-old_w)/dt , dfw )*dx - inner( 0.5*(fp+old_fp) , dfw )*dx
    Res += inner( (p-old_p)/dt , dfp )*dx + inner( 0.5*(fw+old_fw) , dfp )*dx 
    Res += viscosity * inner( grad(fp) , grad(dfp) )*dx
    Res += inner(fw, dw )*dx + inner(fp, dp )*dx - derivative(H, q, dq)
    
    q.assign(old_q)
    solve(Res == 0, q, [bc1,bc2])
    
    # Compute energies
    E_kin = assemble(e_kinetic*dx)
    E_pot = assemble(e_potential*dx)
    return q,E_kin,E_pot
    return q

# Evolve solution t -> t + dt
t = 0
energies = []
old_q = interpolate(initial, W)
for i in range(n_steps+1):
    t+=dt
    q,E_kin,E_pot = evolve(old_q, dt)    
    old_q.assign(q)
    # ... here one would do all the postprocessing, plotting etc
    output(q,t)
    print("iter: %d t: %4.1f E_kin: %8.6f E_pot: %8.6f" % (i,t,E_kin,E_pot))
    energies.append([i,t,E_kin,E_pot,E_kin+E_pot])

E = list(zip(*energies))
t = E[1]
plt.plot(t, E[2], label="E_kin")
plt.plot(t, E[3], label="E_pot")
plt.plot(t, E[4], label="E_tot")
plt.legend(loc="upper right")
plt.xlabel("time")
plt.ylabel("energy")
plt.savefig(output_folder + "energy.png",dpi=300)
