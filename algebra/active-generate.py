from sympy import *
from sympy.printing import print_ccode

# Declare useful symbols
I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar = \
    symbols("I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar")

a, b, a_f, b_f, a_s, b_s, a_fs, b_fs, kappa, beta = \
symbols("a_pt, b_pt, a_f_pt, b_f_pt, a_s_pt, b_s_pt, a_fs_pt, b_fs_pt, kappa_pt, beta_pt")

lambda_a, f_XB, I4_f_e = symbols("lambda_a, f_XB, I_4_f_e")

# Strain energy function in terms of the invariants of the right
# Cauchy-Green tensor for the passive tissue
psi_iso =  a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1)
psi_vol = kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1))

# Elastic energy stored in the cross-bridges
psi_XB = Rational(2, 3)*f_XB*(I4_f_e**Rational(3, 2) - Rational(3, 2)*I4_f_e + Rational(1, 2))

# Identity Matrix
I = eye(3)

# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])
n0 = Matrix([0, 0, 1])

# Compute the elastic components of the stress given the deformation
# gradient
def background_elasticity(C):

    # Modified right Cauchy-Green tensor
    C_bar = J**(-Rational(2, 3))*(C)

    # Define the second Piola-Kirchhoff stress in terms of the invariants
    S_bar =   2*(diff(psi_iso, I1_bar) + diff(psi_iso, I2_bar))*I \
            - 2*diff(psi_iso, I2_bar)*C_bar \
            + 2*diff(psi_iso, I4_f_bar)*(f0*f0.T) \
            + 2*diff(psi_iso, I4_s_bar)*(s0*s0.T) \
            + diff(psi_iso, I8_fs_bar)*(f0*s0.T + s0*f0.T) \
            + diff(psi_iso, I8_fn_bar)*(f0*n0.T + n0*f0.T)
    S_bar_contract_C = sum([sum([S_bar[a, b]*C[a, b]
                                 for a in range(3)]) for b in range(3)])
    Dev_S_bar = S_bar - Rational(1, 3)*(S_bar_contract_C)*C.inv()

    S_iso_inf = J**(-Rational(2, 3))*Matrix(Dev_S_bar)
    S_vol_inf = J*diff(psi_vol, J)*C.inv()

    # Substitute the current values of the invariants
    substitutions = {I8_fn_bar: (f0.T*C_bar*n0)[0]}
    print_ccode(C_bar.trace(), assign_to="real I1_bar")
    print_ccode(Rational(1, 2)*(I1_bar*I1_bar - (C_bar*C_bar).trace()), assign_to="real I2_bar")
    print_ccode((f0.T*C_bar*f0)[0], assign_to="real I4_f_bar")
    print_ccode((s0.T*C_bar*s0)[0], assign_to="real I4_s_bar")
    print_ccode((f0.T*C_bar*s0)[0], assign_to="real I8_fs_bar")
    print_ccode((f0.T*C_bar*n0)[0], assign_to="real I8_fn_bar")

    return S_vol_inf, S_iso_inf


def crossbridge_elasticity(C, lambda_a):

    F_a = Matrix([[lambda_a, 0, 0], [0, 1/sqrt(lambda_a), 0], [0, 0, 1/sqrt(lambda_a)]])
    C_e = F_a.inv().T*C*F_a.inv()
    S_XB = 2*F_a.inv()*diff(psi_XB, I4_f_e)*(f0*f0.T)*F_a.inv().T

    # Substitute the current values of the invariants
    print_ccode((f0.T*C_e*f0)[0], assign_to="real I4_f_e")

    return S_XB

c11, c12, c13, c21, c22, c23, c31, c32, c33 = \
symbols("cloc11, cloc12, cloc13, cloc21, cloc22, cloc23, cloc31, cloc32, cloc33")

C = Matrix([[c11, c12, c13],
            [c21, c22, c23],
            [c31, c32, c33]])

print "// Invariants for background material"
S_vol_inf, S_iso_inf = background_elasticity(C)

print
print "// Invariant for crossbridge elasticity"
S_XB = crossbridge_elasticity(C, lambda_a)

print
print "// Isochoric part of the background elasticity"
for i in range(3):
    for j in range(3):
        print_ccode(S_iso_inf[i, j], assign_to="S_loc_iso(%i, %i)" % (i + 1, j + 1))

print
print "// Volumetric part of the background elasticity"
for i in range(3):
    for j in range(3):
        print_ccode(S_vol_inf[i, j], assign_to="S_loc_vol(%i, %i)" % (i + 1, j + 1))

print
print "// Cross-bridge elasticity"
for i in range(3):
    for j in range(3):
        print_ccode(S_XB[i, j], assign_to="S_XB(%i, %i)" % (i + 1, j + 1))
