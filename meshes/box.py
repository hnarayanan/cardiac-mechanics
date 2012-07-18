from dolfin import *

mesh = Box(-3, 3, -3, 3, -3, 3, 21, 21, 21)
file = File("box.pvd")
file << mesh
