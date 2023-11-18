# compare with src/zoo/testproblems/ROL_Rosenbrock.hpp

import numpy as np
from pyrol import Objective
from pyrol.vectors import NumPyVector

from .TestProblem import TestProblem


class ObjectiveRosenbrock(Objective):
   
    def __init__(self, alpha = 100.0):
        self.alpha = alpha
        super().__init__()

    def value(self, x, tol):
        v = np.sum(self.alpha*(x[::2]**2 - x[1::2])**2 + (x[::2] - 1)**2)
        return v
    
    def gradient(self, g, x, tol):
        g[::2]  =  4*self.alpha*(x[::2]**2 - x[1::2])*x[::2] + 2*(x[::2] - 1)
        g[1::2] = -2*self.alpha*(x[::2]**2 - x[1::2])

    def hessVec(self, hv, v, x, tol):
        h11 =  4*self.alpha*(3*x[::2]**2 - x[1::2]) + 2
        h12 = -4*self.alpha*x[::2]
        h22 =  2*self.alpha
        hv[::2]  = h11*v[::2] + h12*v[1::2]
        hv[1::2] = h12*v[::2] + h22*v[1::2]

    def invHessVec(self, hv, v, x, tol):
        h11 =  4*self.alpha*(3*x[::2]**2 - x[1::2]) + 2
        h12 = -4*self.alpha*x[::2]
        h22 =  2*self.alpha
        hv[::2]  = (+h22*v[::2] - h12*v[1::2])/(h11*h22-h12*h12)
        hv[1::2] = (-h12*v[::2] + h11*v[1::2])/(h11*h22-h12*h12)


class Rosenbrock(TestProblem):

    def __init__(self):
        self.n = 100
   
    def getObjective(self):
        return ObjectiveRosenbrock()
    
    def getInitialGuess(self):
        xnp = np.empty(self.n, dtype='float64')
        xnp[::2]  = -1.2
        xnp[1::2] =  1.0
        return NumPyVector(xnp)

    def getSolution(self, i = 0):
        xnp = np.ones(self.n, dtype='float64')
        return NumPyVector(xnp)