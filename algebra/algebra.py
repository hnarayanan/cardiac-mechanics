from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
I, B, C, gamma, f0, s0, n0 = symbols("I, B, C, gamma, f0, s0, n0")
sigma, p = symbols("sigma, p")

# Material parameters for Figure 7
a    =  0.500 #kPa
b    =  8.023
a_f  = 18.472 #kPa
b_f  = 16.026
a_s  =  2.481 #kPa
b_s  = 11.120
a_fs =  0.356 #kPa
b_fs = 11.436
p    =  0     #kPa

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
F_fs = Matrix([[1, 0, 0],
               [gamma, 1, 0],
               [0, 0, 1]])
sigma_fs = sigma(F_fs)

F_fn = Matrix([[1, 0, 0],
               [0, 1, 0],
               [gamma, 0, 1]])
sigma_fn = sigma(F_fn)

F_sf = F_fs.T
sigma_sf = sigma(F_sf)

F_sn = Matrix([[1, 0, 0],
               [0, 1, 0],
               [0, gamma, 1]])
sigma_sn = sigma(F_sn)

F_nf = F_fn.T
sigma_nf = sigma(F_nf)

F_ns = F_sn.T
sigma_ns = sigma(F_ns)

gammas = np.arange(0, 0.51, 0.01)
sigma_fs_plot = []
sigma_fn_plot = []
sigma_sf_plot = []
sigma_sn_plot = []
sigma_nf_plot = []
sigma_ns_plot = []

for _gamma in gammas:
    sigma_fs_plot.append(sigma_fs[1].subs({gamma:_gamma}))
    sigma_fn_plot.append(sigma_fn[6].subs({gamma:_gamma}))
    sigma_sf_plot.append(sigma_sf[3].subs({gamma:_gamma}))
    sigma_sn_plot.append(sigma_sn[5].subs({gamma:_gamma}))
    sigma_nf_plot.append(sigma_nf[2].subs({gamma:_gamma}))
    sigma_ns_plot.append(sigma_ns[7].subs({gamma:_gamma}))

plt.plot(gammas, sigma_fs_plot, label='fs')
plt.plot(gammas, sigma_fn_plot, label='fn')
plt.plot(gammas, sigma_sf_plot, label='sf')
plt.plot(gammas, sigma_sn_plot, label='sn')
plt.plot(gammas, sigma_nf_plot, label='nf')
plt.plot(gammas, sigma_ns_plot, label='ns')
plt.legend(loc=(0.75, 0.59))
plt.xlabel("Amount of shear")
plt.ylabel("Shear stress (kPa)")
plt.show()
