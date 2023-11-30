# compare with src/zoo/testproblems/ROL_HS42.hpp

import numpy as np
from pyrol import Objective, Constraint, LinearOperator, LinearConstraint
from pyrol.vectors import NumPyVector

from .TestProblem import TestProblem


class ObjectiveHS42(Objective):

    def value(self, x, tol):
        v = (x[0] - 1)**2 + (x[1] - 2)**2  + (x[2] - 3)**2 + (x[3] - 4)**2
        return v  

    def gradient(self, g, x, tol):
        g[0] = 2*(x[0] - 1)
        g[1] = 2*(x[1] - 2)
        g[2] = 2*(x[2] - 3)
        g[3] = 2*(x[3] - 4)

    def hessVec(self, hv, v, x, tol):
        hv[0] = 2*v[0]
        hv[1] = 2*v[1]
        hv[2] = 2*v[2]
        hv[3] = 2*v[3]


class LinearEqualityOperatorHS42(LinearOperator):

    def apply(self, hv, v, tol):
        hv[0] = v[0]
        

class EqualityConstraintHS42(Constraint):

    def value(self, c, x, tol):
        c[0] = x[2]**2 + x[3]**2 - 2

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = 2*x[2]*v[2] + 2*x[3]*v[3]

    def applyAdjointJacobian(self, ajv, v, x, tol):
        ajv[0] = 0
        ajv[1] = 0
        ajv[2] = 2*x[2]*v[0]
        ajv[3] = 2*x[3]*v[0]

    def applyAdjointHessian(self, ahuv, u, v, x, tol):
        ahuv[0] = 0
        ahuv[1] = 0
        ahuv[2] = 2*v[2]*u[0]
        ahuv[3] = 2*v[3]*u[0]


class HS42(TestProblem):

    def __init__(self):
        self.n = 4

    def getObjective(self):
        return ObjectiveHS42()
    
    def getInitialGuess(self):
        xnp = np.ones(self.n, dtype='float64')
        return NumPyVector(xnp)
    
    def getSolution(self, i = 0):
        xnp = np.array([2, 2, 0.6*(2**0.5), 0.8*(2**0.5)], dtype='float64')
        return NumPyVector(xnp)
    
    def getLinearEqualityConstraint(self):
        A = LinearEqualityOperatorHS42()
        b = NumPyVector(np.full(1, -2, dtype='float64'))
        return LinearConstraint(A, b)
        
    def getLinearEqualityMultiplier(self):
        lnp = NumPyVector(np.zeros(1, dtype='float64'))
        return lnp

    def getEqualityConstraint(self):
        return EqualityConstraintHS42()

    def getEqualityMultiplier(self):
        lnp = NumPyVector(np.zeros(1, dtype='float64'))
        return lnp