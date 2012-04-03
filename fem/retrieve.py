from dolfin import *
from numpy import arange

# Dimensions and mesh density of the domain
width = 1
depth = 1
height = 1
n = 10
mesh = Box(0, width, 0, depth, 0, height, n*width, n*depth, n*height)

# Constants related to time stepping
T = 10.0
dt = T/100
times = arange(dt, T + dt, dt)

visco_series = TimeSeries("../output/visco")
V = VectorFunctionSpace(mesh, "Lagrange", 1)
u = Function(V)

for t_n in times:
    visco_series.retrieve(u.vector(), t_n)
    plot(u, mode='displacement')

interactive()
