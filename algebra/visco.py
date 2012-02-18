from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar = \
    symbols("I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar")
gamma = symbols("gamma")

# Material parameters for Figure 7 in HolzapfelOgden2009
a    =  0.500 #kPa
b    =  8.023
a_f  = 18.472 #kPa
b_f  = 16.026
a_s  =  2.481 #kPa
b_s  = 11.120
a_fs =  0.356 #kPa
b_fs = 11.436

# Material parameters for compressibility
kappa = 2.0e6 #kPa
beta = 9.0

# Strain energy function in terms of the invariants of the right
# Cauchy-Green tensor
psi_iso =  a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1)
psi_vol = kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1))


# Compute the elastic components of the stress given the deformation
# gradient
def elastic_stresses(F):
    # Identity Matrix
    I = Matrix([[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]])

    # Right Cauchy-Green tensor
    C = F.T*F

    # Modified right Cauchy-Green tensor
    C_bar = (F.det())**(-Rational(2, 3))*(C)

    # Reference fibre, sheet and sheet-normal directions
    f0 = Matrix([1, 0, 0])
    s0 = Matrix([0, 1, 0])
    n0 = Matrix([0, 0, 1])

    # Define the second Piola-Kirchhoff stress in terms of the invariants
    S_bar =   2*(diff(psi_iso, I1_bar) + diff(psi_iso, I2_bar))*I - 2*diff(psi_iso, I2_bar)*C_bar \
            + 2*diff(psi_iso, I4_f_bar)*(f0*f0.T) + 2*diff(psi_iso, I4_s_bar)*(s0*s0.T) \
            + diff(psi_iso, I8_fs_bar)*(f0*s0.T + s0*f0.T) \
            + diff(psi_iso, I8_fn_bar)*(f0*n0.T + n0*f0.T)
    S_bar_contract_C = sum([sum([S_bar[a, b]*C[a, b]
                                 for a in range(3)]) for b in range(3)])
    Dev_S_bar = S_bar - Rational(1, 3)*(S_bar_contract_C)*C.inv()

    S_iso_inf = (F.det()**(-Rational(2, 3)))*Matrix(Dev_S_bar)
    S_vol_inf = J*diff(psi_vol, J)*C.inv()

    # Substitute the current values of the invariants
    substitutions = {I1_bar: C_bar.trace(),
                     I2_bar: Rational(1, 2)*(I1_bar*I1_bar - (C_bar*C_bar).trace()),
                     J: F.det(),
                     I4_f_bar: (f0.T*C_bar*f0)[0],
                     I4_s_bar: (s0.T*C_bar*s0)[0],
                     I8_fs_bar: (f0.T*C_bar*s0)[0],
                     I8_fn_bar: (f0.T*C_bar*n0)[0]}
    S_vol_inf = S_vol_inf.subs(substitutions)
    S_iso_inf = S_iso_inf.subs(substitutions)

    return S_vol_inf, S_iso_inf


# Subject the body to a known strain protocol and record the stresses
T = 10.0
dt = 0.1
gamma_T = 0.5
times = np.arange(0, T + dt, dt)

stresses = []
strains = []

# Deformation gradient for a simple shear
F = Matrix([[1, 0, 0],
            [gamma, 1, 0],
            [0, 0, 1]])

S_vol_inf, S_iso_inf = elastic_stresses(F)
S = S_vol_inf + S_iso_inf

S = F*S*F.T

for t in times:
    _gamma = t/T*gamma_T
    strains.append(_gamma)
    stresses.append(S.subs({gamma:_gamma})[1, 0])

plt.plot(strains, stresses)
plt.xlabel("Shear strain")
plt.ylabel("Shear stress (kPa)")
plt.show()
