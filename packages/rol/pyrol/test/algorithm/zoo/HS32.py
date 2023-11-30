# compare with src/zoo/testproblems/ROL_HS32.hpp

import numpy as np
from pyrol import Objective, Constraint, LinearOperator, LinearConstraint, Bounds
from pyrol.vectors import NumPyVector

from .TestProblem import TestProblem


class ObjectiveHS32(Objective):

    def value(self, x, tol):
        v = (x[0] + 3*x[1] + x[2])**2 + 4*(x[0] - x[1])**2
        return v  

    def gradient(self, g, x, tol):
        g[0] = 2*(x[0] + 3*x[1] + x[2]) + 8*(x[0] - x[1])
        g[1] = 6*(x[0] + 3*x[1] + x[2]) - 8*(x[0] - x[1])
        g[2] = 2*(x[0] + 3*x[1] + x[2])

    def hessVec(self, hv, v, x, tol):
        hv[0] = 10*v[0] -  2*v[1] + 2*v[2]
        hv[1] = -2*v[0] + 26*v[1] + 6*v[2]
        hv[2] =  2*v[0] +  6*v[1] + 2*v[2]


class LinearEqualityOperatorHS32(LinearOperator):

    def apply(self, hv, v, tol):
        hv[0] = -v[0] - v[1] - v[2]

    def applyAdjoint(self, hv, v, tol):
        hv[0] = -v[0]
        hv[1] = -v[0]
        hv[2] = -v[0]
        

class InequalityConstraintHS32(Constraint):

    def value(self, c, x, tol):
        c[0] = 6*x[1] + 4*x[2] - x[0]**3 - 3

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = - 3*x[0]**2*v[0] + 6*v[1] + 4*v[2]

    def applyAdjointJacobian(self, ajv, v, x, tol):
        ajv[0] = -3*x[0]**2*v[0]
        ajv[1] =  6*v[0]
        ajv[2] =  4*v[0]

    def applyAdjointHessian(self, ahuv, u, v, x, tol):
        ahuv[0] = -6*u[0]*x[0]*v[0]
        ahuv[1] =  0
        ahuv[2] =  0


class HS32(TestProblem):

    def __init__(self):
        self.n = 3

    def getObjective(self):
        return ObjectiveHS32()
    
    def getInitialGuess(self):
        xnp = np.array([0.1, 0.7, 0.2], dtype='float64')
        return NumPyVector(xnp)
    
    def getSolution(self, i = 0):
        xnp = np.array([0, 0, 1], dtype='float64')
        return NumPyVector(xnp)
    
    def getLinearEqualityConstraint(self):
        A = LinearEqualityOperatorHS32()
        b = NumPyVector(np.ones(1, dtype='float64'))
        return LinearConstraint(A, b)
        
    def getLinearEqualityMultiplier(self):
        lnp = NumPyVector(np.zeros(1, dtype='float64'))
        return lnp

    def getInequalityConstraint(self):
        return InequalityConstraintHS32()

    def getInequalityMultiplier(self):
        lnp = NumPyVector(np.zeros(1, dtype='float64'))
        return lnp
    
    def getBoundConstraint(self):
        lower = NumPyVector(np.zeros(self.n, dtype='float64'))
        return Bounds(lower, True)