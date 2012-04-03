# Copyright (C) 2012 Harish Narayanan

# Library imports and settings
from dolfin import *
from numpy import array, arange
parameters["form_compiler"]["name"] = "sfc"

# Dimensions and mesh density of the domain
width = 1
depth = 1
height = 1
n = 10
mesh = Box(0, width, 0, depth, 0, height, n*width, n*depth, n*height)

# Reference fibre, sheet and sheet-normal directions
f0 = Constant((1, 0, 0))
s0 = Constant((0, 1, 0))
n0 = Constant((0, 0, 1))

# Material parameters for Figure 7 in HolzapfelOgden2009
a    =  Constant(0.500)   #kPa
b    =  Constant(8.023)
a_f  =  Constant(16.472)  #kPa
b_f  =  Constant(16.026)
a_s  =  Constant(2.481)   #kPa
b_s  =  Constant(11.120)
a_fs =  Constant(0.356)   #kPa
b_fs =  Constant(11.436)

# Material parameters for compressibility
kappa = Constant(2.0e3)   #kPa
beta  = Constant(9.0)

# Strain energy functions for the passive myocardium
def psi_iso_inf(I1_bar, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar):
   return(a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1))

def psi_vol_inf(J):
    return(kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1)))

# Define kinematic measures in terms of the displacement
def kinematics(u):
    I = Identity(u.cell().d)    # Identity tensor
    F = I + grad(u)             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor
    J = variable(det(F))        # Jacobian
    C_bar = J**(-2.0/3.0)*C     # Modified right Cauchy-Green tensor

    # Principle isotropic invariants
    I1_bar = variable(tr(C_bar))
    I2_bar = variable(0.5*(tr(C_bar)**2 - tr(C_bar*C_bar)))

    # Anisotropic (quasi) invariants
    I4_f_bar = variable(inner(f0, C_bar*f0))
    I4_s_bar = variable(inner(s0, C_bar*s0))
    I8_fs_bar = variable(inner(f0, C_bar*s0))
    I8_fn_bar = variable(inner(f0, C_bar*n0))

    return [I, F, C, J, C_bar, I1_bar, I2_bar, \
            I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar]

# Define the elastic response of the material
# Isochoric part of the second Piola-Kirchhoff stress
def S_iso_inf(u):
    # Define useful kinematic measures
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)

    # Strain energy functions
    psi_iso = psi_iso_inf(I1_bar, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar)

    # Define the second Piola-Kirchhoff stress in terms of the invariants
    S_bar =   2*(diff(psi_iso, I1_bar) + diff(psi_iso, I2_bar))*I \
            - 2*diff(psi_iso, I2_bar)*C_bar \
            + 2*diff(psi_iso, I4_f_bar)*outer(f0, f0) \
            + 2*diff(psi_iso, I4_s_bar)*outer(s0, s0) \
            + diff(psi_iso, I8_fs_bar)*(outer(f0, s0) + outer(s0, f0)) \
            + diff(psi_iso, I8_fn_bar)*(outer(f0, n0) + outer(n0, f0))
    Dev_S_bar = S_bar - (1.0/3.0)*inner(S_bar, C)*inv(C)
    S_iso_inf = J**(-2.0/3.0)*Dev_S_bar
    return(S_iso_inf)

# Volumetric part of the second Piola-Kirchhoff stress
def S_vol_inf(u):
    # Define useful kinematic measures
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)
    psi_vol = psi_vol_inf(J)
    S_vol_inf = J*diff(psi_vol, J)*inv(C)
    return(S_vol_inf)

# First Piola-Kirchhoff stress
def P(u):
    # Define useful kinematic measures
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)
    return (F*(S_iso_inf(u) + S_vol_inf(u)))

# Cauchy stress
def sigma(u):
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)
    return(1/J*P(u)*F.T)

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
shear = Expression(("0.0", "gamma*depth", "0.0"), gamma=0.0, depth=depth)
hold_back = DirichletBC(V, hold, back)
shear_front = DirichletBC(V, shear, front)
bcs = [hold_back, shear_front]

# Define the problem in variational form
F = inner(P(u), grad(v))*dx
J = derivative(F, u, du)

# Constants related to time stepping
T = 10.0
dt = T/100
gamma_max = 0.5
times = arange(dt, T + dt, dt)

#displacement_file = File("../output/displacement.pvd")
#stress_file = File("../output/stress.pvd")
visco_series = TimeSeries("../output/visco")
applied_gamma = 0.0

for t_n in times:

    gamma_n = gamma_max*sin(2*t_n/T*float(pi))
    shear.gamma = gamma_n
    solve(F == 0, u, bcs, J=J)
    visco_series.store(u.vector(), t_n)

#    displacement_file << u
#    stress = project(sigma(u)[0][1], Q) #fs
    # stress = project(sigma(u)[0][2], Q) #fn
    # stress = project(sigma(u)[1][0], Q) #sf
    # stress = project(sigma(u)[1][2], Q) #sn
    # stress = project(sigma(u)[2][0], Q) #nf
    # stress = project(sigma(u)[2][1], Q) #ns
#   center = (depth/2.0, width/2.0, height/2.0)
#    print "stress-strain:", applied_gamma, stress(center)
#    stress_file << stress







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
# shear = Expression(("gamma*height", "0.0", "0.0"), gamma=0.0, height=height)
# hold_bottom = DirichletBC(V, hold, bottom)
# shear_top = DirichletBC(V, shear, top)
# bcs = [hold_bottom, shear_top]

# nf
# shear = Expression(("0.0", "gamma*height", "0.0"), gamma=0.0, height=height)
# hold_bottom = DirichletBC(V, hold, bottom)
# shear_top = DirichletBC(V, shear, top)
# bcs = [hold_bottom, shear_top]


# ffc_options = {
#     "quadrature_degree": 5,
#     "eliminate_zeros": True,
#     "precompute_basis_const": True,
#     "precompute_ip_const": True
#     # "optimize": True
# }
