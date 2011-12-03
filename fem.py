# Copyright (C) 2011 Harish Narayanan

# Library imports and settings
from dolfin import *
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
a    =  0.059 #kPa
b    =  8.023
a_f  = 18.472 #kPa
b_f  = 16.026
a_s  =  2.481 #kPa
b_s  = 11.120
a_fs =  0.216 #kPa
b_fs = 11.436
p    =  a     #kPa

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

# Function spaces
V = VectorFunctionSpace(mesh, "Lagrange", 1)
Q = FunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
#B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
#T  = Constant((0.1, 0.0, 0.0))  # Traction force on the boundary

# Boundary conditions
left_condition   = "x[0] == 0.0 && on_boundary"
right_condition  = "x[0] == %g && on_boundary" % width
back_condition   = "x[1] == 0.0 && on_boundary"
front_condition  = "x[1] == %g && on_boundary" % depth
bottom_condition = "x[2] == 0.0 && on_boundary"
top_condition    = "x[2] == %g && on_boundary" % height

left, right = compile_subdomains([left_condition, right_condition])
back, front = compile_subdomains([back_condition, front_condition])
bottom, top = compile_subdomains([bottom_condition, top_condition])


hold = Expression(("0.0"))
pull = Expression(("strain*height"), height=height, strain=0.0)

hold_left   = DirichletBC(V.sub(0), hold, left)
hold_back   = DirichletBC(V.sub(1), hold, back)
hold_bottom = DirichletBC(V.sub(2), hold, bottom)
#pull_right    = DirichletBC(V.sub(0), pull, right)
#pull_front    = DirichletBC(V.sub(1), pull, front)
pull_top    = DirichletBC(V.sub(2), pull, top)


bcs = [hold_bottom, pull_top, hold_left, hold_back]

F = inner(sigma(u), grad(v))*dx
J = derivative(F, u, du)

displacement_file = File("output/displacement.pvd")
stress_file = File("output/stress.pvd")
applied_strain = 0.0

while applied_strain <= 2.0:
    pull.strain = applied_strain
    solve(F == 0, u, bcs, J=J,
          form_compiler_parameters=ffc_options)
    applied_strain = applied_strain + 0.01
    displacement_file << u
    stress = project(sigma(u)[2][2], Q)
#    stress = project(sigma(u)[0][0], Q)
#    stress = project(sigma(u)[1][1], Q)
#    center = (width/2.0, depth/2.0, height/2.0)
#    print applied_strain
    stress_file << stress




#    plot(u, mode = "displacement")
# Plot and hold solution
#plot(u, mode = "displacement")

#interactive()

