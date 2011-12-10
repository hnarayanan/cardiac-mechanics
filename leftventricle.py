from dolfin import *

mesh = Mesh("meshes/leftventricle.xml")
facet_regions = MeshFunction("uint", mesh,
                             "meshes/leftventricle_facet_region.xml")
regions = {"Endocardium": 1, "Base": 2, "Epicardium": 3}

V = VectorFunctionSpace(mesh, "CG", 1)

hold = Constant((0.0, 0.0, 0.0))
bcs = DirichletBC(V, hold, facet_regions, regions["Base"])

u = TrialFunction(V)
v = TestFunction(V)
f = Constant((0.0, 0.0, 0.0))
p = Constant(1.0)
n = FacetNormal(mesh)

E  = 10.0
nu = 0.3

mu    = E / (2.0*(1.0 + nu))
lmbda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))

def sigma(v):
    return 2.0*mu*sym(grad(v)) + \
        lmbda*tr(sym(grad(v)))*Identity(v.cell().d)


dss = ds[facet_regions]
a = inner(sigma(u), grad(v))*dx
L = inner(f, v)*dx + inner(p*n, v)*dss(regions["Endocardium"])

u = Function(V)
solve(a == L, u, bcs)
plot(u, mode="displacement", interactive=True)
