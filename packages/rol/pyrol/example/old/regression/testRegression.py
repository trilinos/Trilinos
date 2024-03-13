import numpy as np
from numpy.random import default_rng

from pyrol import *
from definitions import myVector, ForwardProblem, Loss


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
