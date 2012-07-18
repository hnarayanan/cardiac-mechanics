from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

l0 = Constant((1, 0, 0))
t0 = Constant((0, 1, 0))

mesh = UnitCube(3, 3, 3)
V = VectorFunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(u, v)*dx

pos = Expression(("x[0]", "x[1]", "x[2]"))

beta = 0.3
v_min = -80
v_max = +40
dv = 2

# Save solution in VTK format
file = File("displacement.pvd");


def gamma_l(v):
    return -beta*(v - v_min)/(v_max - v_min + v)

def gamma_t(v):
    return -gamma_l(v)/(1 + gamma_l(v))

voltages = np.arange(v_min, v_max, dv)

gamma_l_store = []

for v_n in voltages:
    gamma_l_store.append(gamma_t(v_n))

# plt.plot(voltages, gamma_l_store)
# plt.show()

    I = Identity(3)
    F_a = I + gamma_l(v_n)*outer(l0, l0) + gamma_t(v_n)*outer(t0, t0)

    L = inner(F_a*pos - pos, v)*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u)

    file << u;

    # Plot solution
    #plot(u, interactive=True, mode='displacement')
