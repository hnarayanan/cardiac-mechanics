# Copyright (C) 2012 Harish Narayanan

# Library imports and settings
from dolfin import *
from numpy import array, arange

parameters["form_compiler"]["name"] = 'sfc'
ffc_options = {
    "quadrature_degree": 5,
    "eliminate_zeros": True,
    "precompute_basis_const": True,
    "precompute_ip_const": True
    # "optimize": True
}

# Material parameters for Figure 7 in HolzapfelOgden2009
a    =  Constant(0.059)   #kPa
b    =  Constant(8.023)
a_f  =  Constant(18.472)  #kPa
b_f  =  Constant(16.026)
a_s  =  Constant(2.481)   #kPa
b_s  =  Constant(11.120)
a_fs =  Constant(0.216)   #kPa
b_fs =  Constant(11.436)

# Material parameters for compressibility
kappa = Constant(2.0e3)   #kPa
beta  = Constant(9.0)

# Parameters related to time-stepping
T = 1.0
dt = T/10
gamma_max = 0.5

# Strain energy functions for the passive myocardium
# Isochoric part
def psi_iso_inf(I1_bar, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar):
    return(a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1))

# Volumetric part
def psi_vol_inf(J):
    return(kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1)))

# Reference fibre, sheet and sheet-normal directions
f0 = Constant((1, 0, 0))
s0 = Constant((0, 1, 0))
n0 = Constant((0, 0, 1))

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

# Define the elastic response of the material in terms of the second
# Piola-Kirchhoff stress tensor
def S(u):
    # Define useful kinematic measures
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)

    # Strain energy functions
    psi_iso = psi_iso_inf(I1_bar, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar) # Isochoric part
    psi_vol = psi_vol_inf(J)                                                # Volumetric part

    # Define the isochoric part of the second Piola-Kirchhoff stress
    # in terms of the invariants
    S_bar =   2*(diff(psi_iso, I1_bar) + diff(psi_iso, I2_bar))*I \
            - 2*diff(psi_iso, I2_bar)*C_bar \
            + 2*diff(psi_iso, I4_f_bar)*outer(f0, f0) \
            + 2*diff(psi_iso, I4_s_bar)*outer(s0, s0) \
            + diff(psi_iso, I8_fs_bar)*(outer(f0, s0) + outer(s0, f0)) \
            + diff(psi_iso, I8_fn_bar)*(outer(f0, n0) + outer(n0, f0))
    Dev_S_bar = S_bar - (1.0/3.0)*inner(S_bar, C)*inv(C)
    S_iso_inf = J**(-2.0/3.0)*Dev_S_bar

    # Define the volumetric part of the second Piola-Kirchhoff stress
    # in terms of the Jacobian
    S_vol_inf = J*diff(psi_vol, J)*inv(C)

    # Return the total second Piola-Kirchhoff stress
    return(S_iso_inf + S_vol_inf)


# Cauchy stress
def sigma(u):
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)
    return(1/J*P(u)*F.T)

# Dimensions and mesh density of the domain
width = 1
depth = 1
height = 1
n = 5
mesh = Box(0, width, 0, depth, 0, height, n*width, n*depth, n*height)

# Function spaces
scalar = FunctionSpace(mesh, "Lagrange", 1)
vector = VectorFunctionSpace(mesh, "Lagrange", 1)
tensor = TensorFunctionSpace(mesh, "Lagrange", 1)

# Functions
du = TrialFunction(vector)            # Incremental displacement
v  = TestFunction(vector)             # Test function
u  = Function(vector)                 # Displacement from previous iteration

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

# Simple shear along the fs plane
shear = Expression(("0.0", "gamma*depth", "0.0"), gamma=0.0, depth=depth)
hold_back = DirichletBC(vector, hold, back)
shear_front = DirichletBC(vector, shear, front)
bcs = [hold_back, shear_front]


# Define the time range
times = arange(dt, T + dt, dt)

# Subject the body to a known strain protocol and record the stresses
for t_n in times:

    print "I was here at %g s" % t_n

    # Compute current shear strain and update the boundary condition
    gamma_n = gamma_max*t_n/T
    shear.gamma = gamma_n

    # Define the variational form for the problem
    F = inner((Identity(u.cell().d) + grad(u))*S(u), grad(v))*dx
    J = derivative(F, u, du)

    # Solve the boundary value problem
    solve(F == 0, u, bcs, J=J)

#    stress = project(S_n(u)[0][1], scalar) #fs
#    center = (depth/2.0, width/2.0, height/2.0)
#    print "stress-strain:", gamma_n, stress(center)
