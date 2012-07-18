from sympy import *

# Declare useful symbols
q_1, q_2, E_1, E_2, f_1, f_2, lmbda_opt, xi, v = \
    symbols("q_1, q_2, E_1, E_2, f_1, f_2, lmbda_opt, xi, v")
lmbda, lmbda_a, lmbda_c, alpha_1, alpha_2, alpha_3, alpha_4 = \
    symbols("lambda, lambda_a, lambda_c, alpha(1), alpha(2), alpha(3), alpha(4)")

psi_1 = q_1/q_2*(exp(q_1*(lmbda - 1)) - 1 + q_2*(1 - lmbda))
psi_2 = (alpha_3*E_1 + alpha_4*E_2)*Rational(1, 2)*(lmbda_c - 1)**2
N = 1#exp(-(lmbda_a - lmbda_opt)**2/(2*xi**2))
C = (f_1*alpha_3 + f_2*alpha_4)*N
P_a = -f_1*alpha_3*v*N

P = diff(psi_1, lmbda) + N/lmbda_a*diff(psi_2, lmbda_c)
print P

lmbda_a_dot = 1/C*(P_a - diff(N, lmbda_a)*psi_2 + (lmbda_c/lmbda_a)*N*diff(psi_2, lmbda_c))
print lmbda_a_dot

# f0 = Matrix([1, 0, 0])
# s0 = Matrix([0, 1, 0])
# n0 = Matrix([0, 0, 1])

# F_a = lmbda_a*f0*f0.T + 1/sqrt(lmbda_a)*(s0*s0.T + n0*n0.T)
# dF_dlmbda_a = f0*f0.T - Rational(1, 2)*lmbda_a**Rational(3, 2)*(s0*s0.T + n0*n0.T)

# psi_2 = Rational(2, 3)*(alpha_3*E_1 + alpha_4*E_2)*()

# lmbd_a_dot_3d = 1/C*(P_a - diff(N, lmbda_a)*psi_2 + 2*F_a.inv().T*F_a.inv()*[]*F_a.inv().T:dF_dlmbda_a

#    S_bar_contract_C = sum([sum([S_bar[a, b]*C[a, b]
 #                                for a in range(3)]) for b in range(3)])
