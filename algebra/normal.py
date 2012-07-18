from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
I, B, C, gamma, f0, s0, n0 = symbols("I, B, C, gamma, f0, s0, n0")
sigma, p = symbols("sigma, p")

# Material parameters for Figure 7
a    =  1.0 #kPa
b    =  1.0
a_f  =  0.0 #kPa
b_f  =  0.0
a_s  =  0.0 #kPa
b_s  =  0.0
a_fs =  0.0 #kPa
b_fs =  0.0
p    =  a   #kPa

# a    =  0.500 #kPa
# b    =  8.023
# a_f  = 18.472 #kPa
# b_f  = 16.026
# a_s  =  2.481 #kPa
# b_s  = 11.120
# a_fs =  0.356 #kPa
# b_fs = 11.436
# p    =  0     #kPa

def sigma(F):
    # Cauchy-Green tensors
    B = F*F.T
    C = F.T*F

    # Principle isotropic invariants
    I1 = C.trace()
    I2 = Rational(1, 2)*(I1*I1 - (C*C).trace())
    I3 = C.det()

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

    # Cauchy stress
    sigma =   a*exp(b*(I1 - 3))*B - p*I \
        + 2*a_f*(I4_f - 1)*exp(b_f*(I4_f - 1)**2)*(f*f.T) \
        + 2*a_s*(I4_s - 1)*exp(b_s*(I4_s - 1)**2)*(s*s.T) \
        + a_fs*I8_fs*exp(b_fs*I8_fs**2)*(f*s.T + s*f.T)

    return(sigma)

    # Returning the Piola-Kirchhoff tensor instead results
    # in something strange
    # return(I3*sigma*F.inv().T)

# Identity Matrix
I = Matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

# Deformation gradients and corresponding stresses
F_ff = Matrix([[1 + gamma, 0, 0],
               [0, 1/sqrt(1 + gamma), 0],
               [0, 0, 1/sqrt(1 + gamma)]])
sigma_ff = sigma(F_ff)

gammas = np.arange(0, 0.11, 0.01)
sigma_ff_plot = []

for _gamma in gammas:
    sigma_ff_plot.append(sigma_ff[0, 0].subs({gamma:_gamma}))

plt.plot(gammas, sigma_ff_plot, label='ff')
plt.legend(loc=(0.75, 0.59))
plt.xlabel("Shear strain")
plt.ylabel("Shear stress (kPa)")
plt.show()
