import numpy as np
from numpy.random import default_rng

from PyROL import *
from PyROL.vectors import npVector


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


def main():
    
    a, b = 1, 0  # the "true" parameter values from which we generate data
    n = 100      # the size of our data set
    
    # Configure parameter list.  ################
    params = getParametersFromXmlFile("input.xml")
    
    # Set the output stream.  ###################
    stream = getCout()
    
    # Create predictors.  #######################
    rng = default_rng(23040901)
    X = rng.uniform(-1, 1, n)
    
    # Create the forward problem.  ##############
    constraint = ForwardProblem(X)

    # Set up vectors.  ##########################
    u = myVector(np.empty(n))
    z = myVector(np.empty(2))
    g = z.dual()  # gradient with respect to z
    l = u.dual()  # adjoint variables
    
    # Create responses.  ########################
    Y = u.clone()
    constraint.solve(u.clone(), Y, myVector(np.array([a, b])), np.nan)
    noise = 0.1*rng.normal(0, 1, n)
    Y.values += noise
    
    # Set up the problem.  ######################
    objective = Loss(Y)
    reducedObjective = Reduced_Objective_SimOpt(objective, constraint, u, z, l)
    problem = Problem(reducedObjective, z, g)
    problem.check(True, stream)
    
    # Solve.  ###################################
    solver = Solver(problem, params)
    solver.solve(stream)
    
    # Check the solution.  ######################
    target, _, _, _ = np.linalg.lstsq(constraint.X, objective.Y.values, rcond=None)
    np.testing.assert_allclose(z.values, target)
    print('Test passed: Correctly optimized!')


if __name__ == "__main__":
    main()