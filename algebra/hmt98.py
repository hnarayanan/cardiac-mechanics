from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

beta_0 = 1.45
beta_1 = 1.95
beta_2 = 0.31
T_ref = 125
n_ref = 4.25
pC_50_ref = 5.33

def T_0(lmbda, Ca):

    pC_50 = pC_50_ref*(1 + beta_2*(lmbda - 1))
    C_50 = 10**(6 - pC_50)
    n = n_ref*(1 + beta_1*(lmbda - 1))
    z_ss = Ca**n/(Ca**n + C_50**n)
    z = z_ss

    T_0 = T_ref*(1 + beta_0*(lmbda - 1))*z
    return(T_0)


concentrations = np.arange(0.01, 100 + 0.01, 0.01)
lmbdas = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3]

for lmbda in lmbdas:
    T_store = []
    for Ca in concentrations:
        T_store.append(T_0(lmbda, Ca)/T_ref)
    plt.semilogx(concentrations, T_store, label='$\lambda$ = %.2f' % lmbda)

plt.legend(loc='best')
plt.show()

