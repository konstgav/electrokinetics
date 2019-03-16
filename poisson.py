# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square

def SolvePoisson():
    ngsglobals.msg_level = 1

    # generate a triangular mesh of mesh-size 0.2
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.02))

    # H1-conforming finite element space
    fes = H1(mesh, order=3, dirichlet=[1,2,3,4])

    # define trial- and test-functions
    u = fes.TrialFunction()
    v = fes.TestFunction()

    # the right hand side
    f = LinearForm(fes)
    f += SymbolicLFI(v)

    # the bilinear-form 
    a = BilinearForm(fes, symmetric=True)
    a += SymbolicBFI(grad(u)*grad(v))

    a.Assemble()
    f.Assemble()

    # the solution field 
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # plot the solution (netgen-gui only)
    Draw (gfu)

    exact = 16*x*(1-x)*y*(1-y)
    print ("L2-error:", sqrt (Integrate ( (gfu-exact)*(gfu-exact), mesh)))

    return gfu, mesh