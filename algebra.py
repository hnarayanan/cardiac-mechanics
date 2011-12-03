from sympy import *
import numpy as np
import matplotlib.pyplot as plt


print("-"*72)
print("Passive response of the myocardium")
print("-"*72)

# Declare useful symbols
I, B, C, gamma, f0, s0, n0 = symbols("I, B, C, gamma, f0, s0, n0")
a, b, a_f, b_f, a_s, b_s, a_fs, b_fs = \
    symbols("a, b, a_f, b_f, a_s, b_s, a_fs, b_fs")
sigma, p = symbols("sigma, p")

# Identity Matrix
I = Matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

# Deformation gradient
F = Matrix([[1, gamma, 0],
            [0, 1, gamma],
            [0, 0, 1]])

# Cauchy-Green tensors
B = F*F.T
C = F.T*F

# Principle isotropic invariants
I1 = B.trace()
I2 = Rational(1, 2)*(I1*I1 - (B*B).trace())
I3 = B.det()

# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])

# Anisotropic (quasi) invariants
I4_f = (f0.T*C*f0)[0]
I4_s = (s0.T*C*s0)[0]
I8_fs = (f0.T*C*s0)[0]

# Current fibre, sheet and sheet-normal directions
f = F*f0
s = F*s0

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

# Cauchy stress
sigma =   a*exp(b*(I1 - 3))*B - p*I \
        + 2*a_f*(I4_f - 1)*exp(b_f*(I4_f - 1)**2)*(f*f.T) \
        + 2*a_s*(I4_s - 1)*exp(b_s*(I4_s - 1)**2)*(s*s.T) \
        + a_fs*I8_fs*exp(b_fs*I8_fs**2)*(f*s.T + s*f.T)

gammas = np.arange(0, 0.51, 0.01)
sigma0 = []
for _gamma in gammas:
    sigma0.append(sigma[0].subs({gamma:_gamma}))

plt.plot(gammas, sigma0)
plt.show()
