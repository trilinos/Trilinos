from PyROL.vectors.npVector import npVector as vector_type
from PyROL import *
from PyROL.PyROL.ROL import getParametersFromXmlFile
import numpy as np


class norm2Obj(Objective):
    def __init__(self, H, g, c=0):
        self.H = H
        self.g = g
        self.c = c
        super().__init__()

    def value(self, x, tol):
        tmp = x.clone()
        x.zero()
        self.H.apply(tmp, x, tol)
        tmp.scale(0.5)
        tmp.plus(self.g)
        return x.apply(tmp) + self.c

    def gradient(self, g, x, tol):
        self.H.apply(g, x, tol)
        g.plus(self.g)

    def hessVec(self, hv, v, x, tol):
        self.H.apply(hv, v, tol)

    def invHessVec(self, hv, v, x, tol):
        self.H.applyInverse(hv, v, tol)


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
obj = norm2Obj(op, g)
c = Constraint()
a = vector_type.full(10, 1.)
b = vector_type.full(10, 1.)
a.scale(2.)
np.testing.assert_equal(a.norm(), 2*np.sqrt(10))
a.zero()
np.testing.assert_equal(a.norm(), 0.)

a.axpy(1., b)
np.testing.assert_equal(a.norm(), np.sqrt(10))
np.testing.assert_equal(a.dot(b), 10.)

b.setScalar(2.)
a = b
np.testing.assert_equal(a.apply(b), 40.)
np.testing.assert_equal(b.norm(), 2*np.sqrt(10))
np.testing.assert_equal(b[0], 2.)
b[0:2] = [1., 3.]
np.testing.assert_equal(b[0:2], [1., 3.])

print(obj.value(a, 1e-8))

op.apply(b,a,1e-8)
print(b[0:3])
g = vector_type.full(10, 0.)

params = getParametersFromXmlFile("input.xml")

problem = Problem(obj, x)
#solver = Solver(problem, params)

#print(params)

obj_q = QuadraticObjective(op, g)
step = TrustRegionStep(params)
status = StatusTest(params)
algo = Algorithm(step,status,False)

print("before optimization: "+str(obj_q.value(x, 1e-8)))
algo.run_void(x=x, obj=obj_q)
print("after optimization: "+str(obj_q.value(x, 1e-8)))
print(g[:])
print(x[:])
