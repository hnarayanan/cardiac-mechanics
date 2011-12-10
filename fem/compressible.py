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
mesh = Box(0, 0, 0, depth, width, height, n*depth, n*width, n*height)

# Material parameters
c1 = Constant(1.0)
c2 = Constant(0.1)
kappa = Constant(1.0e6)

def P(u):
    I = Identity(u.cell().d)    # Identity tensor
    F = I + grad(u)             # Deformation gradient
    J = variable(det(F))        # Jacobian
    C = F.T*F                   # Right Cauchy-Green tensor

    # Decoupling the volumetric and isochoric parts
    F_bar = J**(-1.0/3.0)*F
    C_bar = J**(-2.0/3.0)*F
    C_bar_inv = inv(C_bar)

    # Principle isotropic invariants
    I1_bar = variable(tr(C_bar))
    I2_bar = variable(0.5*(tr(C_bar)**2 - tr(C_bar*C_bar)))
    I3_bar = det(C_bar) # = 1

    # Define the material model
    psi_iso = c1*(I1_bar - 3) + c2*(I1_bar - 3)**2
    psi_vol = kappa*(0.25*(J**2 - 1 - 2*ln(J)))

    gamma1_bar = 2*(diff(psi_iso, I1_bar) + I1_bar*diff(psi_iso, I2_bar))
    gamma2_bar = -2*diff(psi_iso, I2_bar)
    S_bar = 2*(gamma1_bar*I + gamma2_bar*C_bar)
    S_iso = J**(-2.0/3.0)*as_matrix([[-0.333333333333333*C_bar[0][1]*C_bar_inv[0][0]*S_bar[0][1] - 0.333333333333333*C_bar[0][2]*C_bar_inv[0][0]*S_bar[0][2] - 0.333333333333333*C_bar[1][0]*C_bar_inv[0][0]*S_bar[1][0] - 0.333333333333333*C_bar[1][1]*C_bar_inv[0][0]*S_bar[1][1] - 0.333333333333333*C_bar[1][2]*C_bar_inv[0][0]*S_bar[1][2] - 0.333333333333333*C_bar[2][0]*C_bar_inv[0][0]*S_bar[2][0] - 0.333333333333333*C_bar[2][1]*C_bar_inv[0][0]*S_bar[2][1] - 0.333333333333333*C_bar[2][2]*C_bar_inv[0][0]*S_bar[2][2] + S_bar[0][0]*(-0.333333333333333*C_bar[0][0]*C_bar_inv[0][0] + 1.0), -0.333333333333333*C_bar[0][0]*C_bar_inv[1][0]*S_bar[0][0] - 0.333333333333333*C_bar[0][2]*C_bar_inv[1][0]*S_bar[0][2] - 0.333333333333333*C_bar[1][1]*C_bar_inv[1][0]*S_bar[1][1] - 0.333333333333333*C_bar[1][2]*C_bar_inv[1][0]*S_bar[1][2] - 0.333333333333333*C_bar[2][0]*C_bar_inv[1][0]*S_bar[2][0] - 0.333333333333333*C_bar[2][1]*C_bar_inv[1][0]*S_bar[2][1] - 0.333333333333333*C_bar[2][2]*C_bar_inv[1][0]*S_bar[2][2] + S_bar[0][1]*(-0.333333333333333*C_bar[0][1]*C_bar_inv[1][0] + 0.5) + S_bar[1][0]*(-0.333333333333333*C_bar[1][0]*C_bar_inv[1][0] + 0.5), -0.333333333333333*C_bar[0][0]*C_bar_inv[2][0]*S_bar[0][0] - 0.333333333333333*C_bar[0][1]*C_bar_inv[2][0]*S_bar[0][1] - 0.333333333333333*C_bar[1][0]*C_bar_inv[2][0]*S_bar[1][0] - 0.333333333333333*C_bar[1][1]*C_bar_inv[2][0]*S_bar[1][1] - 0.333333333333333*C_bar[1][2]*C_bar_inv[2][0]*S_bar[1][2] - 0.333333333333333*C_bar[2][1]*C_bar_inv[2][0]*S_bar[2][1] - 0.333333333333333*C_bar[2][2]*C_bar_inv[2][0]*S_bar[2][2] + S_bar[0][2]*(-0.333333333333333*C_bar[0][2]*C_bar_inv[2][0] + 0.5) + S_bar[2][0]*(-0.333333333333333*C_bar[2][0]*C_bar_inv[2][0] + 0.5)], [-0.333333333333333*C_bar[0][0]*C_bar_inv[0][1]*S_bar[0][0] - 0.333333333333333*C_bar[0][2]*C_bar_inv[0][1]*S_bar[0][2] - 0.333333333333333*C_bar[1][1]*C_bar_inv[0][1]*S_bar[1][1] - 0.333333333333333*C_bar[1][2]*C_bar_inv[0][1]*S_bar[1][2] - 0.333333333333333*C_bar[2][0]*C_bar_inv[0][1]*S_bar[2][0] - 0.333333333333333*C_bar[2][1]*C_bar_inv[0][1]*S_bar[2][1] - 0.333333333333333*C_bar[2][2]*C_bar_inv[0][1]*S_bar[2][2] + S_bar[0][1]*(-0.333333333333333*C_bar[0][1]*C_bar_inv[0][1] + 0.5) + S_bar[1][0]*(-0.333333333333333*C_bar[1][0]*C_bar_inv[0][1] + 0.5), -0.333333333333333*C_bar[0][0]*C_bar_inv[1][1]*S_bar[0][0] - 0.333333333333333*C_bar[0][1]*C_bar_inv[1][1]*S_bar[0][1] - 0.333333333333333*C_bar[0][2]*C_bar_inv[1][1]*S_bar[0][2] - 0.333333333333333*C_bar[1][0]*C_bar_inv[1][1]*S_bar[1][0] - 0.333333333333333*C_bar[1][2]*C_bar_inv[1][1]*S_bar[1][2] - 0.333333333333333*C_bar[2][0]*C_bar_inv[1][1]*S_bar[2][0] - 0.333333333333333*C_bar[2][1]*C_bar_inv[1][1]*S_bar[2][1] - 0.333333333333333*C_bar[2][2]*C_bar_inv[1][1]*S_bar[2][2] + S_bar[1][1]*(-0.333333333333333*C_bar[1][1]*C_bar_inv[1][1] + 1.0), -0.333333333333333*C_bar[0][0]*C_bar_inv[2][1]*S_bar[0][0] - 0.333333333333333*C_bar[0][1]*C_bar_inv[2][1]*S_bar[0][1] - 0.333333333333333*C_bar[0][2]*C_bar_inv[2][1]*S_bar[0][2] - 0.333333333333333*C_bar[1][0]*C_bar_inv[2][1]*S_bar[1][0] - 0.333333333333333*C_bar[1][1]*C_bar_inv[2][1]*S_bar[1][1] - 0.333333333333333*C_bar[2][0]*C_bar_inv[2][1]*S_bar[2][0] - 0.333333333333333*C_bar[2][2]*C_bar_inv[2][1]*S_bar[2][2] + S_bar[1][2]*(-0.333333333333333*C_bar[1][2]*C_bar_inv[2][1] + 0.5) + S_bar[2][1]*(-0.333333333333333*C_bar[2][1]*C_bar_inv[2][1] + 0.5)], [-0.333333333333333*C_bar[0][0]*C_bar_inv[0][2]*S_bar[0][0] - 0.333333333333333*C_bar[0][1]*C_bar_inv[0][2]*S_bar[0][1] - 0.333333333333333*C_bar[1][0]*C_bar_inv[0][2]*S_bar[1][0] - 0.333333333333333*C_bar[1][1]*C_bar_inv[0][2]*S_bar[1][1] - 0.333333333333333*C_bar[1][2]*C_bar_inv[0][2]*S_bar[1][2] - 0.333333333333333*C_bar[2][1]*C_bar_inv[0][2]*S_bar[2][1] - 0.333333333333333*C_bar[2][2]*C_bar_inv[0][2]*S_bar[2][2] + S_bar[0][2]*(-0.333333333333333*C_bar[0][2]*C_bar_inv[0][2] + 0.5) + S_bar[2][0]*(-0.333333333333333*C_bar[2][0]*C_bar_inv[0][2] + 0.5), -0.333333333333333*C_bar[0][0]*C_bar_inv[1][2]*S_bar[0][0] - 0.333333333333333*C_bar[0][1]*C_bar_inv[1][2]*S_bar[0][1] - 0.333333333333333*C_bar[0][2]*C_bar_inv[1][2]*S_bar[0][2] - 0.333333333333333*C_bar[1][0]*C_bar_inv[1][2]*S_bar[1][0] - 0.333333333333333*C_bar[1][1]*C_bar_inv[1][2]*S_bar[1][1] - 0.333333333333333*C_bar[2][0]*C_bar_inv[1][2]*S_bar[2][0] - 0.333333333333333*C_bar[2][2]*C_bar_inv[1][2]*S_bar[2][2] + S_bar[1][2]*(-0.333333333333333*C_bar[1][2]*C_bar_inv[1][2] + 0.5) + S_bar[2][1]*(-0.333333333333333*C_bar[2][1]*C_bar_inv[1][2] + 0.5), -0.333333333333333*C_bar[0][0]*C_bar_inv[2][2]*S_bar[0][0] - 0.333333333333333*C_bar[0][1]*C_bar_inv[2][2]*S_bar[0][1] - 0.333333333333333*C_bar[0][2]*C_bar_inv[2][2]*S_bar[0][2] - 0.333333333333333*C_bar[1][0]*C_bar_inv[2][2]*S_bar[1][0] - 0.333333333333333*C_bar[1][1]*C_bar_inv[2][2]*S_bar[1][1] - 0.333333333333333*C_bar[1][2]*C_bar_inv[2][2]*S_bar[1][2] - 0.333333333333333*C_bar[2][0]*C_bar_inv[2][2]*S_bar[2][0] - 0.333333333333333*C_bar[2][1]*C_bar_inv[2][2]*S_bar[2][1] + S_bar[2][2]*(-0.333333333333333*C_bar[2][2]*C_bar_inv[2][2] + 1.0)]])
    p = diff(psi_vol, J)
    S_vol = J*p*inv(C)
    S = S_iso + S_vol
    return(F*S)

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

# nf
shear = Expression(("0.0", "gamma*height", "0.0"), gamma=0.0, height=height)
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
    stress = project(P(u)[2][1], Q) #ns
    print "stress-strain:", applied_gamma, max(stress.vector().array())
    stress_file << stress
