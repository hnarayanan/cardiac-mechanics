from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
# Local symbols
I1, I2, I3, I4_f, I4_s, I8_fs, I8_fn = symbols("I1, I2, I3, I4_f, I4_s, I8_fs, I8_fn")
gamma = symbols("gamma")

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

# Identity Matrix
I = Matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

# Define the strain energy function
def strain_energy():
    psi  =   a/(2*b)*exp(b*(I1 - 3)) \
           + a_f/(2*b_f)*(exp(b_f*(I4_f - 1)**2) - 1) \
           + a_s/(2*b_s)*(exp(b_s*(I4_s - 1)**2) - 1) \
           + a_fs/(2*b_fs)*(exp(b_fs*I8_fs**2) - 1)
    return(psi)

def sigma(F):

    # Cauchy-Green tensors
    B = F*F.T
    C = F.T*F

    # Reference fibre, sheet and sheet-normal directions
    f0 = Matrix([1, 0, 0])
    s0 = Matrix([0, 1, 0])
    n0 = Matrix([0, 0, 1])

    # Current fibre, sheet and sheet-normal directions
    f = F*f0
    s = F*s0
    n = F*n0

    # Define the Cauchy stress
    psi = strain_energy()
    sigma = - p*I + 2*diff(psi, I1)*B \
            + 2*diff(psi, I2)*(I1*B - B**2) \
            + 2*diff(psi, I4_f)*(f*f.T) \
            + 2*diff(psi, I4_s)*(s*s.T) \
            + diff(psi, I8_fs)*(f*s.T + s*f.T) \
            + diff(psi, I8_fn)*(f*n.T + n*f.T)

    # Principle isotropic invariants
    _I1 = C.trace()
    _I2 = Rational(1, 2)*(I1*I1 - (C*C).trace())
    _I3 = C.det()

    # Anisotropic (quasi) invariants
    _I4_f = (f0.T*C*f0)[0]
    _I4_s = (s0.T*C*s0)[0]
    _I8_fs = (f0.T*C*s0)[0]
    _I8_fn = (f0.T*C*n0)[0]


    return(sigma.subs({I1: _I1, I2:_I2, I3:_I3,
                       I4_f:_I4_f, I4_s:_I4_s,
                       I8_fs: _I8_fs, I8_fn:_I8_fn}))

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
plt.xlabel("Shear strain")
plt.ylabel("Shear stress (kPa)")
plt.show()
