from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Declare useful symbols
t, s, tau, alpha, eps_0 = symbols('t s tau alpha eps_0')
E, E_inf, E_0 = symbols('E, E_inf, E_0')

# alpha = 1/tau*integrate(exp(-(t - s)/tau)*eps(s), (s, -oo, t))
alpha = 1/tau*integrate(exp(-(t - s)/tau)*eps_0, (s, 0, t))
E = E_0 - E_inf
sigma = E_0*eps_0 - E*alpha
#tau = -sigma.subs({t: 0})/diff(sigma, t).subs({t: 0})*(E/E_0)

times = np.arange(0, 10.05, 0.05)
sigma_plot = []

for _time in times:
    sigma_plot.append(sigma.subs({E_inf: 1.0, E_0: 1.1, \
                                  eps_0: 2, tau: 2.50,  \
                                  t: _time}))

plt.plot(times, sigma_plot, label='sigma')
plt.show()
