# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry

def SolvePoisson(width, height):
    ngsglobals.msg_level = 1

    # generate a triangular mesh of mesh-size 0.2
    #mesh = Mesh(unit_square.GenerateMesh(maxh=0.02))
    geo = SplineGeometry()
    p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0,0), (width,0), (width, height), (0, height)] ]
    geo.Append (["line", p1, p2],bc=1)
    geo.Append (["line", p2, p3],bc=2)
    geo.Append (["line", p3, p4],bc=3)
    geo.Append (["line", p4, p1],bc=4)
    mesh = Mesh(geo.GenerateMesh(maxh=0.02))
    
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

    return gfu, mesh