#import netgen.gui
from netgen.geom2d import unit_square
from ngsolve import *
from numpy import empty
import matplotlib.pyplot as plt
from matplotlib import cm
from poisson import SolvePoisson
from math import cosh, sin, pi

#TODO: save to file, plot
def GetSolutionOnRectGrid(field, mesh, Nx, Ny):
    res = empty([Ny,Nx],float)
    hx = 1.0/(Nx-1)
    hy = 1.0/(Ny-1)
    for i in range(Nx):
        for j in range(Ny):
            res[j,i] = field(mesh(i*hx, j*hy))
    return res

#TODO: save text data
def PlotAndSave(data, title):
    fig, ax = plt.subplots()
    cax = ax.imshow(data, cmap=cm.coolwarm, origin='lower', extent=(0,1,0,1))
    ax.set_title(title)
    cbar = fig.colorbar(cax)
    plt.savefig('fig.png')
    plt.show()

def GetAnalyticalSolution(Nx, Ny):
    res = empty([Ny,Nx],float)
    hx = 1.0/(Nx-1)
    hy = 1.0/(Ny-1)
    height = 1.
    width = 1.
    pressure = 1.
    for i in range(Nx):
        for j in range(Ny):
            res[j,i] = 0.
            for n in range(1,20,2):
                res[j,i] += (1-cosh(n*pi*(i*hx-width/2.))/cosh(n*pi*width/2./height))*sin(n*pi*j*hy/height)/n**3
            res[j,i] *= 4.*height**2*pressure/(pi**3)
    return res

Nx = 100
Ny = 100
gfu, mesh = SolvePoisson()
numericSolution = GetSolutionOnRectGrid(gfu, mesh, Nx, Ny)
analytSolution = GetAnalyticalSolution(Nx, Ny)
residual = 0.
for i in range(Nx):
    for j in range(Ny):
        residual = max(abs(numericSolution[i,j]-analytSolution[i,j]), residual)
print('Max  residual = {}'.format(residual))

#PlotAndSave(numericSolution, 'Poisson equation solution')
#PlotAndSave(analytSolution, 'Velocity analytical solution')