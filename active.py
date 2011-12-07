# Copyright (C) 2011 Harish Narayanan

# Library imports and settings
from dolfin import *
from numpy import array

parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

# Dimensions and mesh density of the domain
width = 1
depth = 10
height = 1
n = 4
mesh = Box(0, 0, 0, depth, width, height, n*depth, n*width, n*height)
plot(mesh)
interactive()

# Material parameters
# Figure 7
a    = Constant(0.059)  #kPa
b    = Constant(8.023)
a_f  = Constant(18.472) #kPa
b_f  = Constant(16.026)
a_s  = Constant(2.481)  #kPa
b_s  = Constant(11.120)
a_fs = Constant(0.216)  #kPa
b_fs = Constant(11.436)
p    = Constant(0.059)  #kPa

# Reference fibre, sheet and sheet-normal directions
f0 = Constant((1, 0, 0))
s0 = Constant((0, 1, 0))

# Define the material
def sigma(u, gamma):
    I = Identity(u.cell().d)      # Identity tensor
    F = I + grad(u)               # Deformation gradient

    # F_a = I - gamma*outer(f0, f0) # Active strain
    F_e = F*(I + (gamma/(1 - gamma))*outer(f0, f0))

    C_e = F_e.T*F_e               # Right Cauchy-Green tensor
    B_e = F_e*F_e.T               # Left Cauchy-Green tensor

    # Principle isotropic invariants
    I1 = tr(B_e)
    I2 = 0.5*(tr(B_e)**2 - tr(B_e*B_e))
    I3 = det(B_e)

    # Anisotropic (quasi) invariants
    I4_f = inner(f0, C_e*f0)
    I4_s = inner(s0, C_e*s0)
    I8_fs = inner(f0, C_e*s0)

    # Current fibre, sheet and sheet-normal directions
    f = F*f0
    s = F*s0

    # Cauchy stress
    return(  a*exp(b*(I1 - 3))*B_e - p*I \
           + 2*a_f*(I4_f - 1)*exp(b_f*(I4_f - 1)**2)*outer(f, f) \
           + 2*a_s*(I4_s - 1)*exp(b_s*(I4_s - 1)**2)*outer(s, s) \
           + a_fs*I8_fs*exp(b_fs*I8_fs**2)*(outer(f, s) + outer(s, f)))

def P(u, gamma):
    I = Identity(u.cell().d)
    F = I + grad(u)
    return(det(F)*sigma(u, gamma)*inv(F).T)

# Function spaces
V = VectorFunctionSpace(mesh, "Lagrange", 1)
Q = FunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration

# Boundary conditions
back_condition   = "x[0] == 0.0 && on_boundary"
front_condition  = "x[0] == %g && on_boundary" % depth
left_condition   = "x[1] == 0.0 && on_boundary"
right_condition  = "x[1] == %g && on_boundary" % width
bottom_condition = "x[2] == 0.0 && on_boundary"
top_condition    = "x[2] == %g && on_boundary" % height

back, front = compile_subdomains([back_condition, front_condition])
left, right = compile_subdomains([left_condition, right_condition])
bottom, top = compile_subdomains([bottom_condition, top_condition])

hold = Expression(("0.0", "0.0", "0.0"))
hold_back = DirichletBC(V, hold, back)
#hold_front = DirichletBC(V, hold, front)
bcs = [hold_back]

displacement_file = File("output/displacement.pvd")
stress_file = File("output/stress.pvd")
t = 0.0

while t <= 2*DOLFIN_PI:
    active_strain = 0.1*(1 - sin(t - 3*DOLFIN_PI/2))
    gamma = Constant(active_strain)
    F = inner(P(u, gamma), grad(v))*dx
    J = derivative(F, u, du)
    solve(F == 0, u, bcs, J=J,
          form_compiler_parameters=ffc_options)
    t = t + 0.1
    displacement_file << u
    stress = project(sigma(u, gamma)[0][0], Q) #ff
    print "stress-strain:", active_strain, max(stress.vector().array())
    stress_file << stress
