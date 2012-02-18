from dolfin import *

mesh = Mesh("biventricle.xml")
file = File("biventricle.pvd")
file << mesh
