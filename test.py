#import netgen.gui
from netgen.geom2d import unit_square
from ngsolve import *
from numpy import empty
import matplotlib.pyplot as plt
from matplotlib import cm
from poisson import SolvePoisson
from math import cosh, sin, pi

#TODO: save to file, plot
def GetSolutionOnRectGrid(field, mesh, width, height, Nx, Ny):
    res = empty([Ny,Nx],float)
    hx = width/(Nx-1)
    hy = height/(Ny-1)
    for i in range(Nx):
        for j in range(Ny):
            res[j,i] = field(mesh(i*hx, j*hy))
    return res

#TODO: save text data
def PlotAndSave(data, title, width, height):
    fig, ax = plt.subplots()
    cax = ax.imshow(data, cmap=cm.coolwarm, origin='lower', extent=(0,width,0,height))
    ax.set_title(title)
    cbar = fig.colorbar(cax)
    plt.savefig(title+'.png')
    plt.show()

def GetAnalyticalSolution(width, height, Nx, Ny):
    res = empty([Ny,Nx],float)
    hx = width/(Nx-1)
    hy = height/(Ny-1)
    pressure = 1.
    for i in range(Nx):
        for j in range(Ny):
            res[j,i] = 0.
            for n in range(1,40,2):
                res[j,i] += (1-cosh(n*pi*(i*hx-width/2.)/height)/cosh(n*pi*width/2./height))*sin(n*pi*j*hy/height)/n**3
            res[j,i] *= 4.*height**2*pressure/(pi**3)
    return res

Nx = 100
Ny = 100
width = 1.
height = 2.
gfu, mesh = SolvePoisson(width, height)
numericSolution = GetSolutionOnRectGrid(gfu, mesh, width, height, Nx, Ny)
analytSolution = GetAnalyticalSolution(width, height, Nx, Ny)
residual = 0.
imax = 0
jmax = 0
for i in range(Nx):
    for j in range(Ny):
        res = abs(numericSolution[i,j]-analytSolution[i,j])
        if (res>residual):
            residual = res
            imax = i
            jmax = j
print('Max  residual = {}, posision i = {}, j = {}'.format(residual, imax, jmax))

PlotAndSave(numericSolution, 'Poisson equation solution', width, height)
PlotAndSave(analytSolution, 'Velocity analytical solution', width, height)
