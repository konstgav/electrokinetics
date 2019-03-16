from numpy import e, tanh, log, linspace, sinh, abs, vstack, array, zeros, cosh
from pylab import show, plot, legend, xlabel, ylabel
from scipy.integrate import solve_bvp

def InfinitePlane():
    x = linspace(0,5)
    phi0 = 1.
    phiDebyeHuckle = phi0*e**(-x)
    phi = 2*log((1+e**(-x)*tanh(0.25*phi0))/(1-e**(-x)*tanh(0.25*phi0)))
    plot(x,phiDebyeHuckle,label='Debye-Huckle linearization')
    plot(x,phi, label = 'Poisson-Boltzmann equation solution')
    legend()
    show()

def fun(x, y):
    return vstack((y[1], k**2*sinh(y[0])))

def bc(ya, yb):
    return array([abs(ya[0]-phi0), abs(yb[0]-phi0)])

def SymmetryPlanes():
    x = linspace(-1,1,100)
    
    phiDebyeHuckle = phi0*cosh(k*x)/cosh(k)
    y_a = zeros((2, x.size))

    res_a = solve_bvp(fun, bc, x, y_a)
    y_plot_a = res_a.sol(x)[0]

    plot(x,phiDebyeHuckle,label='Debye-Huckle linearization')
    plot(x,y_plot_a, label = 'Poisson-Boltzmann equation solution')
    xlabel('x')
    ylabel(r'$ varphi $')
    legend()
    show()

phi0 =6.
k = 4.
SymmetryPlanes()