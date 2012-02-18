from sympy import *
from sympy.abc import t

tau = symbols('tau')
Q = Function('Q')
S_iso = Function('S_iso')

eqn = Derivative(Q(t), t) + Q(t)/tau - Derivative(S_iso(t), t)
print(dsolve(eqn, Q(t)))
