from PyROL.vectors.npVector import npVector as vector_type
from PyROL import *
from PyROL.PyROL.ROL import getParametersFromXmlFile
from PyROL.PyROL.ROL import Bounds_double_t, Problem_double_t, Solver_double_t
import numpy as np

# Matrix from rol/example/quadratic/example_01.cpp
class matrix(LinearOperator):
    def __init__(self, dim):
        self.dim = dim
        super().__init__()

    def apply(self, Hv, v, tol):
        for i in range(0, self.dim):
            Hv[i] = 2.*v[i]
            if i > 0:
                Hv[i] -= v[i-1]
            if i < self.dim - 1:
                Hv[i] -= v[i+1]

op = matrix(10)
g = vector_type.full(10, 0.)
x = vector_type.full(10, 1.)

# passes:
# bnd = Bounds_double_t(g.clone(), x.clone())
# fails:
bnd = Bounds_double_t(vector_type.full(10, 0.), x.clone())
obj = QuadraticObjective(op, g)

params = getParametersFromXmlFile("input.xml")
status = StatusTest(params)

problem = Problem_double_t(obj, x)
problem.addBoundConstraint(bnd)
solver = Solver_double_t(problem, params)
solver.solve(status)
