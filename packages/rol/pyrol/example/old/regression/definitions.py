import numpy as np

from pyrol import *
from pyrol.vectors import npVector


class myVector(npVector):
    
    def norm(self):
        return np.sqrt(self.dot(self))
    
    def __sub__(self, other):
        return myVector(self.values - other.values)


class ForwardProblem(Constraint_SimOpt):
    
    def __init__(self, X):
        self.X = np.reshape(X, (X.size, 1))
        self.X = np.hstack([self.X, np.ones_like(self.X)])
        super().__init__()
        
    def value(self, c, u, z, tol):
        c.values[:] = u.values - self.X@z.values
        
    def applyJacobian_1(self, jv, v, u, z, tol):
        jv.values[:] = v.values
    
    def applyJacobian_2(self, jv, v, u, z, tol):
        jv.values[:] = - self.X@v.values
    
    def applyInverseJacobian_1(self, ijv, v, u, z, tol):
        ijv.values[:] = v.values
        
    def applyAdjointJacobian_1(self, ajv, v, u, z, tol):
        ajv.values[:] = v.values
    
    def applyAdjointJacobian_2(self, ajv, v, u, z, tol):
        ajv.values[:] = - self.X.T@v.values
    
    def applyInverseAdjointJacobian_1(self, iajv, v, u, z, tol):
        iajv.values[:] = v.values
        
    def solve(self, c, u, z, tol):
        u.values[:] = self.X@z.values
        self.value(c, u, z, tol)


class Loss(Objective_SimOpt):
    
    def __init__(self, Y):
        self.Y = Y
        super().__init__()
        
    def value(self, u, z, tol):
        return (self.Y - u).norm()**2/2
    
    def gradient_1(self, g, u, z, tol):
        g.values[:] = u.values - self.Y.values
        
    def gradient_2(self, g, u, z, tol):
        g.zero()
        
    def hessVec_11(self, hv, v, u, z, tol):
        hv.values[:] = v.values
        
    def hessVec_12(self, hv, v, u, z, tol):
        hv.zero()
        
    def hessVec_21(self, hv, v, u, z, tol):
        hv.zero()
        
    def hessVec_22(self, hv, v, u, z, tol):
        hv.zero()
