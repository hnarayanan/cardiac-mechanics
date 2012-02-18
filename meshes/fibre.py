from dolfin import *

mesh = Mesh("leftventricle.xml")
fibre = Expression(("0.0", "-b*b/(-a*a)*(x[0]/x[1])*x[0]", "0.0"), a = 2.5, b = 1.1, c = 1.1)
plot(fibre, mesh=mesh, interactive=True)
