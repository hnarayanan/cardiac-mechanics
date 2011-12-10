from dolfin import *

mesh = Mesh("meshes/leftventricle.xml")
facet_regions = MeshFunction("uint", mesh,
                             "meshes/leftventricle_facet_region.xml")
regions = {"Endocardium": 1, "Base": 2, "Epicardium": 3}

V = VectorFunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

hold = Constant((0.0, 0.0, 0.0))
bcs = DirichletBC(V, hold, facet_regions, regions["Base"])

E  = 10.0
nu = 0.3

mu    = E / (2.0*(1.0 + nu))
lmbda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))

def sigma(v):
    return 2.0*mu*sym(grad(v)) + \
        lmbda*tr(sym(grad(v)))*Identity(v.cell().d)


f = Constant((0.0, 0.0, 0.0))
n = FacetNormal(mesh)
dss = ds[facet_regions]

t = 0.0
displacement_file = File("output/displacement.pvd")

a = inner(sigma(u), grad(v))*dx
L1 = inner(f, v)*dx

u = Function(V)

while t <= 2*DOLFIN_PI:

    p = Constant((1 - sin(t - 3*DOLFIN_PI/2) - 0.5))
    L2 = inner(p*n, v)*dss(regions["Endocardium"])
    L = L1 + L2
    solve(a == L, u, bcs)
    t = t + 0.1
    displacement_file << u
