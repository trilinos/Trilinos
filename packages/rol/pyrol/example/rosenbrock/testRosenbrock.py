import numpy as np

from pyrol import *
from pyrol.vectors import npVector as myVector


class RosenbrockObjective(Objective):
    
    def __init__(self):
        self.alpha = 100
        super().__init__()
        
    def value(self, x, tol):
        v = self.alpha*(x.values[0]**2 - x.values[1])**2 + (x.values[0] - 1)**2
        return v
    
    def gradient(self, g, x, tol):
        g[0]  =  2*self.alpha*(x.values[0]**2 - x.values[1])*2*x.values[0] 
        g[0] +=  2*(x.values[0] - 1)
        g[1]  = -2*self.alpha*(x.values[0]**2 - x.values[1])
        
    def hessVec(self, hv, v, x, tol):
        h11 = 12*self.alpha*x.values[0]**2 - 4*self.alpha*x.values[1] + 2
        h12 = -4*self.alpha*x.values[0]
        h22 = 2*self.alpha
        hv[0] = h11*v[0] + h12*v[1]
        hv[1] = h12*v[0] + h22*v[1]


def main():
    
    # Configure parameter list.  ################
    params = getParametersFromXmlFile("input.xml")
    # How do we initialize an instance without the command above?
    # ROL::ParameterList parlist;
    params.sublist("General").sublist("Secant").set("Use as Hessian", False)
    # params.sublist("Step").set("Type", "Trust Region");
    params.sublist("Step").sublist("Trust Region").set("Subproblem Solver", "Truncated CG")
    
    # Set the output stream.  ###################
    stream = getCout()
    
    # Set up vectors.  ##########################
    x = myVector(np.array([-3., -4.]))
    g = x.dual()
    
    # Set up the problem.  ######################
    objective = RosenbrockObjective()
    problem = Problem(objective, x, g)
    problem.check(True, stream)
    
    # Solve.  ###################################
    solver = Solver(problem, params)
    solver.solve(stream)
    
    # Check the solution.  ######################
    target = np.array([1., 1.])
    np.testing.assert_allclose(x.values, target)
    print('Test passed: Correctly optimized!')


if __name__ == "__main__":
    main()
