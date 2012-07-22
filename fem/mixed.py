# Copyright (C) 2012 Harish Narayanan

# Library imports and settings
from dolfin import *
from numpy import array, arange

parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {
    "quadrature_degree": 2,
    "eliminate_zeros": True,
    "precompute_basis_const": True,
    "precompute_ip_const": True
    # "optimize": True
}

# Material constants
a    =  Constant(0.500)   #kPa
b    =  Constant(8.023)
kappa = Constant(2.0e3)   #kPa
beta  = Constant(9.0)

# Strain energy functions
# Isochoric part
def psi_iso(I1_bar):
    return(a/(2*b)*exp(b*(I1_bar - 3)))

# Volumetric part
def psi_vol(J):
    return(kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1)))

# Define the elastic response of the material
# Isochoric part of the second Piola-Kirchhoff stress
def S_iso(I1_bar, J):
    # Define the second Piola-Kirchhoff stress in terms of the invariants
    S_bar =   2*(diff(psi_iso(I1_bar), I1_bar))*I
    Dev_S_bar = S_bar - (1.0/3.0)*inner(S_bar, C)*inv(C)
    S_iso = J**(-2.0/3.0)*Dev_S_bar
    return(S_iso)

# Volumetric part of the second Piola-Kirchhoff stress
def S_vol(J):
    S_vol = J*diff(psi_vol(J), J)*inv(C)
    return(S_vol)

# Define a suitable mesh
n = 10
mesh = UnitCube(n, n, n)

# Function spaces
scalar_DG0 = FunctionSpace(mesh, "DG", 0)
vector_CG1 = VectorFunctionSpace(mesh, "Lagrange", 1)
ME = MixedFunctionSpace([vector_CG1, scalar_DG0, scalar_DG0])

# Functions
dU = TrialFunction(ME)          # Incremental displacement
V  = TestFunction(ME)           # Test function
U  = Function(ME)               # Displacement from previous iteration

u, p, Jbar = split(U)
Jbar = variable(Jbar)
v, q, theta = split(V)
du, dp, dJbar = split(dU)


# Define kinematics
I = Identity(u.cell().d)    # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor
J = variable(det(F))        # Jacobian
C_bar = J**(-2.0/3.0)*C     # Modified right Cauchy-Green tensor
I1_bar = variable(tr(C_bar)) # First invariant

S = S_iso(I1_bar, J)# + S_vol(J)


functional1 = inner(F*S, grad(v))*dx
functional2 = (J - Jbar)*q*dx
functional3 = (diff(psi_vol(Jbar), Jbar) - p)*theta*dx
functional = functional1 + functional2 + functional3

jac = derivative(functional, U, dU)

# Boundary conditions
back_condition   = "x[0] == 0.0 && on_boundary"
front_condition  = "x[0] == 1.0 && on_boundary"

back, front = compile_subdomains([back_condition, front_condition])
hold = Expression(("0.0", "0.0", "0.0"))

shear = Expression(("0.0", "gamma", "0.0"), gamma=0.1)
hold_back = DirichletBC(ME.sub(0), hold, back)
shear_front = DirichletBC(ME.sub(0), shear, front)
bcs = [hold_back, shear_front]

# Solve the boundary value problem
solve(functional == 0, U, bcs, J=jac)

plot(u, interactive=True)
