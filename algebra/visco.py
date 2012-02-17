from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
I1, I2, I3, I4_f, I4_s, I8_fs, I8_fn = symbols("I1, I2, I3, I4_f, I4_s, I8_fs, I8_fn")
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

# Strain energy function in terms of the invariants of the right
# Cauchy-Green tensor
psi = a/(2*b)*exp(b*(I1 - 3)) \
    + a_f/(2*b_f)*(exp(b_f*(I4_f - 1)**2) - 1) \
    + a_s/(2*b_s)*(exp(b_s*(I4_s - 1)**2) - 1) \
    + a_fs/(2*b_fs)*(exp(b_fs*I8_fs**2) - 1)

# Identity Matrix
I = Matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

# Deformation gradient for a simple shear, gamma
F = Matrix([[1, 0, 0],
            [gamma, 1, 0],
            [0, 0, 1]])

# Right Cauchy-Green tensor
C = F.T*F

# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])
n0 = Matrix([0, 0, 1])

# Define the second Piola-Kirchhoff stress in terms of the invariants
S =   2*(diff(psi, I1) + diff(psi, I2))*I - 2*diff(psi, I2)*C \
    + 2*diff(psi, I3)*C.inv() \
    + 2*diff(psi, I4_f)*(f0*f0.T) + 2*diff(psi, I4_s)*(s0*s0.T) \
    + diff(psi, I8_fs)*(f0*s0.T + s0*f0.T) \
    + diff(psi, I8_fn)*(f0*n0.T + n0*f0.T)

# Substitute the current values of the invariants
S = S.subs({I1: C.trace(),
            I2: Rational(1, 2)*(I1*I1 - (C*C).trace()),
            I3: C.det(),
            I4_f: (f0.T*C*f0)[0],
            I4_s: (s0.T*C*s0)[0],
            I8_fs: (f0.T*C*s0)[0],
            I8_fn: (f0.T*C*n0)[0]})

S = F*S*F.T

gammas = np.arange(0, 0.51, 0.01)
S_plot = []

for _gamma in gammas:
    S_plot.append(S[1].subs({gamma:_gamma}))

plt.plot(gammas, S_plot)
plt.xlabel("Shear strain")
plt.ylabel("Shear stress (kPa)")
plt.show()
