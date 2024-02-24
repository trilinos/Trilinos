# compare with src/zoo/testproblems/ROL_HS41.hpp

import numpy as np
from pyrol import Objective, LinearOperator, LinearConstraint, Bounds
from pyrol.vectors import NumPyVector

from .TestProblem import TestProblem


class ObjectiveHS41(Objective):

    def value(self, x, tol):
        v = 2 - x[0]*x[1]*x[2]
        return v  

    def gradient(self, g, x, tol):
        g[0] = -x[1]*x[2]
        g[1] = -x[0]*x[2]
        g[2] = -x[0]*x[1]
        g[3] = 0

    def hessVec(self, hv, v, x, tol):
        hv[0] = -x[2]*v[1] - x[1]*v[2]
        hv[1] = -x[2]*v[0] - x[0]*v[2]
        hv[2] = -x[1]*v[0] - x[0]*v[1]
        hv[3] = 0


class LinearEqualityOperatorHS41(LinearOperator):

    def apply(self, hv, v, tol):
        hv[0] = v[0] + 2*v[1] + 2*v[2] - v[3]
        
    def applyAdjoint(self, hv, v, tol):
        hv[0] =   v[0] 
        hv[1] = 2*v[0]
        hv[2] = 2*v[0]
        hv[3] = - v[0]


class HS41(TestProblem):

    def __init__(self):
        self.n = 4

    def getObjective(self):
        return ObjectiveHS41()
    
    def getInitialGuess(self):
        xnp = np.full(self.n, 2, dtype='float64')
        return NumPyVector(xnp)
    
    def getSolution(self, i = 0):
        xnp = np.array([2/3, 1/3, 1/3, 2], dtype='float64')
        return NumPyVector(xnp)
    
    def getLinearEqualityConstraint(self):
        A = LinearEqualityOperatorHS41()
        b = NumPyVector(np.zeros(1, dtype='float64'))
        return LinearConstraint(A, b)

    def getLinearEqualityMultiplier(self):
        lnp = NumPyVector(np.zeros(1, dtype='float64'))
        return lnp

    def getBoundConstraint(self):
        lower = NumPyVector(np.zeros(self.n, dtype='float64'))
        upper = NumPyVector(np.array([1, 1, 1, 2], dtype='float64'))
        return Bounds(lower, upper)