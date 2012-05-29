from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

# Declare useful symbols
I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar = \
    symbols("I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar")
gamma = symbols("gamma")

a, b, a_f, b_f, a_s, b_s, a_fs, b_fs = \
    symbols('a, b, a_f, b_f, a_s, b_s, a_fs, b_fs')
kappa, beta = symbols('kappa, beta')

material_parameters = {a: 0.500, b: 8.023, a_f: 16.472, b_f: 16.026, \
                       a_s: 2.481, b_s: 11.120, a_fs: 0.356, b_fs: 11.436, \
                       kappa: 2.0e6, beta:9.0}

# Strain energy function in terms of the invariants of the right
# Cauchy-Green tensor
psi_iso =  a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1)
psi_vol = kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1))
# psi_vol = kappa*(J*ln(J) - J + 1)

# Identity Matrix
I = eye(3)

# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])
n0 = Matrix([0, 0, 1])

# Compute the elastic components of the stress given the deformation
# gradient
def elastic_stresses(F):

    # Right Cauchy-Green tensor
    C = F.T*F

    # Modified right Cauchy-Green tensor
    C_bar = (F.det())**(-Rational(2, 3))*(C)

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

    S_iso_inf = (F.det()**(-Rational(2, 3)))*Matrix(Dev_S_bar)
    S_vol_inf = J*diff(psi_vol, J)*C.inv()

    # Substitute the current values of the invariants
    kinematics = {I1_bar: C_bar.trace(),
                     I2_bar: Rational(1, 2)*(I1_bar*I1_bar - (C_bar*C_bar).trace()),
                     J: F.det(),
                     I4_f_bar: (f0.T*C_bar*f0)[0],
                     I4_s_bar: (s0.T*C_bar*s0)[0],
                     I8_fs_bar: (f0.T*C_bar*s0)[0],
                     I8_fn_bar: (f0.T*C_bar*n0)[0]}
    substitutions = dict(kinematics.items() + material_parameters.items())
    S_vol_inf = S_vol_inf.subs(substitutions)
    S_iso_inf = S_iso_inf.subs(substitutions)

    return S_vol_inf, S_iso_inf


# Constants related to time stepping
T = 10.0
dt = T/100
gamma_max = 0.1
times = np.arange(dt, T + dt, dt)

# Constants related to viscoelasticity
tau = 1.5
beta_inf = 2.0
xi = -dt/(2*tau)

# Deformation gradient for tension along the fibre direction
F = Matrix([[1 + gamma, 0, 0],
               [0, 1/sqrt(1 + gamma), 0],
               [0, 0, 1/sqrt(1 + gamma)]])

# Initialize matrices
S_vol_inf_store = [zeros(3)]
S_iso_inf_store = [zeros(3)]
Q_store = [zeros(3)]
S_store = [0.0]
gamma_store = [0.0]

# Analytical values for the elastic stresses in terms of strain
S_vol_inf, S_iso_inf = elastic_stresses(F)

# Subject the body to a known strain protocol and record the stresses
for t_n in times:

    # Load previous elastic stress states
    S_vol_inf_p = S_vol_inf_store[-1]
    S_iso_inf_p = S_iso_inf_store[-1]
    Q_p = Q_store[-1]

    # Compute current strain measures
    # gamma_n = gamma_max*sin(2*t_n/T*float(pi))
    gamma_n = gamma_max*t_n/T
    gamma_store.append(gamma_n)

    # Update stress state
    S_vol_inf_n = S_vol_inf.subs({gamma:gamma_n})
    S_iso_inf_n = S_iso_inf.subs({gamma:gamma_n})
    H_p = exp(xi)*(exp(xi)*Q_p - beta_inf*S_iso_inf_p)
    Q_n = beta_inf*exp(xi)*S_iso_inf_n + H_p
    S_n = S_vol_inf_n + S_iso_inf_n + Q_n

    # Convert to Cauchy stress for comparison with Dokos et al.
    S_n = F.subs({gamma:gamma_n})*S_n*F.T.subs({gamma:gamma_n})

    # Store stress state at current time
    S_vol_inf_store.append(S_vol_inf_n)
    S_iso_inf_store.append(S_iso_inf_n)
    Q_store.append(Q_n)
    S_store.append(S_n[0, 0])

# Plot the current stress-strain curve
plt.plot(gamma_store, S_store, label='ff')


plt.legend(loc=(0.6, 0.6))
plt.xlabel("Strain")
plt.ylabel("Stress (kPa)")
plt.show()
