from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
I1, I2, J, I4_f, I4_s, I8_fs, I8_fn = symbols("I1, I2, J, I4_f, I4_s, I8_fs, I8_fn")
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
psi_iso =  a/(2*b)*exp(b*(I1 - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs**2) - 1)
psi_vol = kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1))

# Identity Matrix
I = Matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

# Deformation gradient for a simple shear, gamma
F = Matrix([[1, 0, 0],
            [gamma, 1, 0],
            [0, 0, 1]])

# Modified right Cauchy-Green tensor
C = (F.det())**(-Rational(2, 3))*(F.T*F)


# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])
n0 = Matrix([0, 0, 1])

# Define the second Piola-Kirchhoff stress in terms of the invariants
S_iso =   2*(diff(psi_iso, I1) + diff(psi_iso, I2))*I - 2*diff(psi_iso, I2)*C \
        + 2*diff(psi_iso, I4_f)*(f0*f0.T) + 2*diff(psi_iso, I4_s)*(s0*s0.T) \
        + diff(psi_iso, I8_fs)*(f0*s0.T + s0*f0.T) \
        + diff(psi_iso, I8_fn)*(f0*n0.T + n0*f0.T)
S_vol = J*diff(psi_vol, J)*C.inv()
S = S_vol + S_iso


# Substitute the current values of the invariants
S = S.subs({I1: C.trace(),
            I2: Rational(1, 2)*(I1*I1 - (C*C).trace()),
            J: F.det(),
            I4_f: (f0.T*C*f0)[0],
            I4_s: (s0.T*C*s0)[0],
            I8_fs: (f0.T*C*s0)[0],
            I8_fn: (f0.T*C*n0)[0]})

gammas = np.arange(0, 0.51, 0.01)
S_plot = []

for _gamma in gammas:
    S_plot.append(S[1].subs({gamma:_gamma}))

plt.plot(gammas, S_plot)
plt.xlabel("Shear strain")
plt.ylabel("Shear stress (kPa)")
plt.show()
