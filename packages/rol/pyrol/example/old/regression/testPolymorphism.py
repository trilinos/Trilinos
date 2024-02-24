import numpy as np
from numpy.random import default_rng

from pyrol import *
from definitions import myVector, ForwardProblem

class Loss(Objective_SimOpt):
    
    def __init__(self, Y):
        self.Y = Y
        super().__init__()
    
    def value_3(self, *args):
        u, z, tol = args
        return (self.Y - u).norm()**2/2

def main():
    
    a, b = 1, 0  # the "true" parameter values from which we generate data
    n = 100      # the size of our data set
        
    # Create predictors.  #######################
    rng = default_rng(23040901)
    X = rng.uniform(-1, 1, n)
    
    # Create the forward problem.  ##############
    constraint = ForwardProblem(X)

    # Set up vectors.  ##########################
    u = myVector(np.empty(n))
    z = myVector(np.empty(2))
    
    # Create responses.  ########################
    Y = u.clone()
    constraint.solve(u.clone(), Y, myVector(np.array([a, b])), np.nan)
    noise = 0.1*rng.normal(0, 1, n)
    Y.values += noise
    
    # Set up the objective.  ####################
    objective = Loss(Y)
    
    # Evaluate the objective.  ##################
    x = Vector_SimOpt(u, z)
    tol = 0.0
    print(objective.value(u, z, tol))
    print(objective.value(x, tol))


if __name__ == "__main__":
    main()
