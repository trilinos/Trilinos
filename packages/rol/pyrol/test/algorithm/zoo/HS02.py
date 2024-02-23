# compare with src/zoo/testproblems/ROL_HS2.hpp

import numpy as np
from pyrol import Objective, LinearOperator, LinearConstraint, Bounds
from pyrol.vectors import NumPyVector

from .Rosenbrock import Rosenbrock


class HS02(Rosenbrock):

    def __init__(self):
        self.n = 2

    def getInitialGuess(self):
        xnp = np.array([-2, 1], dtype='float64')
        return NumPyVector(xnp)
    
    def getSolution(self, i = 0):
        a = (598/1200)**0.5
        b = 400*a**3
        xnp = np.array([2*a*np.cos(np.arccos(1/b)/3), 1.5], dtype='float64')
        return NumPyVector(xnp)

    def getBoundConstraint(self):
        lower = NumPyVector(np.array([-float('inf'), 1.5], dtype='float64'))
        return Bounds(lower, True)