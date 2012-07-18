from dolfin import *
parameters["form_compiler"]["name"] = "ffc"

mesh = UnitSquare(2, 2)
V = VectorFunctionSpace(mesh, "Lagrange", 1)
u  = Function(V)

for j in range(1, 100):
    print j
    term = grad(u)[0][0]
    print term
