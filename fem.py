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
depth = 1
height = 1
n = 5
mesh = Box(0, width, 0, depth, 0, height, n*width, n*depth, n*height)

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
def sigma(u):
    I = Identity(u.cell().d)    # Identity tensor
    F = I + grad(u)             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor
    B = F*F.T                   # Left Cauchy-Green tensor

    # Principle isotropic invariants
    I1 = tr(B)
    I2 = 0.5*(tr(B)**2 - tr(B*B))
    I3 = det(B)

    # Anisotropic (quasi) invariants
    I4_f = inner(f0, C*f0)
    I4_s = inner(s0, C*s0)
    I8_fs = inner(f0, C*s0)

    # Current fibre, sheet and sheet-normal directions
    f = F*f0
    s = F*s0

    # Cauchy stress
    return(  a*exp(b*(I1 - 3))*B - p*I \
           + 2*a_f*(I4_f - 1)*exp(b_f*(I4_f - 1)**2)*outer(f, f) \
           + 2*a_s*(I4_s - 1)*exp(b_s*(I4_s - 1)**2)*outer(s, s) \
           + a_fs*I8_fs*exp(b_fs*I8_fs**2)*(outer(f, s) + outer(s, f)))

def P(u):
    I = Identity(u.cell().d)
    F = I + grad(u)
    return(det(F)*sigma(u)*inv(F).T)

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

# fs
# shear = Expression(("0.0", "gamma*depth", "0.0"), gamma=0.0, depth=depth)
# hold_back = DirichletBC(V, hold, back)
# shear_front = DirichletBC(V, shear, front)
# bcs = [hold_back, shear_front]

# fn
# shear = Expression(("0.0", "0.0", "gamma*depth"), gamma=0.0, depth=depth)
# hold_back = DirichletBC(V, hold, back)
# shear_front = DirichletBC(V, shear, front)
# bcs = [hold_back, shear_front]

# sf
# shear = Expression(("gamma*width", "0.0", "0.0"), gamma=0.0, width=width)
# hold_left = DirichletBC(V, hold, left)
# shear_right = DirichletBC(V, shear, right)
# bcs = [hold_left, shear_right]

# sn
# shear = Expression(("0.0", "0.0", "gamma*width"), gamma=0.0, width=width)
# hold_left = DirichletBC(V, hold, left)
# shear_right = DirichletBC(V, shear, right)
# bcs = [hold_left, shear_right]

# nf
shear = Expression(("gamma*height", "0.0", "0.0"), gamma=0.0, height=height)
hold_bottom = DirichletBC(V, hold, bottom)
shear_top = DirichletBC(V, shear, top)
bcs = [hold_bottom, shear_top]

F = inner(P(u), grad(v))*dx
J = derivative(F, u, du)

displacement_file = File("output/displacement.pvd")
stress_file = File("output/stress.pvd")
applied_gamma = 0.0

while applied_gamma <= 0.50:
    shear.gamma = applied_gamma
    solve(F == 0, u, bcs, J=J,
          form_compiler_parameters=ffc_options)
    applied_gamma = applied_gamma + 0.01
    displacement_file << u
    # stress = project(sigma(u)[0][1], Q) #fs
    # stress = project(sigma(u)[0][2], Q) #fn
    # stress = project(sigma(u)[1][0], Q) #sf
    # stress = project(sigma(u)[1][2], Q) #sn
    stress = project(sigma(u)[2][0], Q) #nf
    print "stress-strain:", applied_gamma, max(stress.vector().array())
    stress_file << stress
