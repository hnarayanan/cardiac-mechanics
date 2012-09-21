from sympy import *
from sympy.printing import print_ccode

# Declare useful symbols
lmbda_a, mu, xi, E_1, E_2, E_3, E_4, XB_PreR, XB_PostR = \
    symbols("lambda_a, mu, xi, E_1, E_2, E_3, E_4, XB_PreR, XB_PostR")
I4_f_e, C, v = symbols("I4_f_e, C, v")

# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])
n0 = Matrix([0, 0, 1])

# Pick a suitable form for the active strain
F_a = lmbda_a*(f0*f0.T) + 1/sqrt(lmbda_a)*(s0*s0.T + n0*n0.T)

# Make other suitable choices for constitutive functions
N = exp(-(lmbda_a - mu)**2/(2*xi**2))
psi_2 = N*Rational(2, 3)*(E_1*XB_PreR + E_2*XB_PostR) \
         *(I4_f_e**Rational(3, 2) - Rational(3, 2) + Rational(1, 2))
P_a = 0 # something
P_a = -N*E_3*XB_PreR*v
C = N*(E_3*XB_PreR + E_4*XB_PostR)

# Define an arbitrary deformation gradient tensor
f11, f12, f13, f21, f22, f23, f31, f32, f33 = \
symbols("f11, f12, f13, f21, f22, f23, f31, f32, f33")

F = Matrix([[f11, f12, f13],
            [f21, f22, f23],
            [f31, f32, f33]])

# Compute the necessary elastic decompositions
F_e = F*F_a.inv()
C_e = F_e.T*F_e


# Compute the more complicted third  term on the right-hand side of the ODE
term_3a = 2*C_e*diff(psi_2, I4_f_e)*(f0*f0.T)*F_a.inv().T
term_3b = Matrix([[diff(F_a[0, 0], lmbda_a), diff(F_a[0, 1], lmbda_a), diff(F_a[0, 2], lmbda_a)], \
                        [diff(F_a[1, 0], lmbda_a), diff(F_a[1, 1], lmbda_a), diff(F_a[1, 2], lmbda_a)], \
                        [diff(F_a[2, 0], lmbda_a), diff(F_a[2, 1], lmbda_a), diff(F_a[2, 2], lmbda_a)]])
term_3  = sum([sum([term_3a[a, b]*term_3b[a, b]
                        for a in range(3)]) for b in range(3)])

# Define the entire right-hand side of the ODE
lmbda_a_dot = 1/C*(P_a - diff(psi_2, lmbda_a) + term_3)

print_ccode((f0.T*C_e*f0)[0], assign_to = "real I4_f_e")
print_ccode(lmbda_a_dot, assign_to = "lambda_a_dot")
