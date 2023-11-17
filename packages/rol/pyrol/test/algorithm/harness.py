import numpy as np

from pyrol import getCout, Solver
from pyrol.vectors import NumPyVector

def harness(mockProblem, parameterList):
    stream = getCout()
    problem, x, xstar = mockProblem.get()
   
    # TO-DO: print problem name

    error = x.clone()
    solver = Solver(problem, parameterList)
    solver.solve(stream)

    e = np.inf
    for xi in xstar:
        error[:] = x[:] - xi[:]
        e = min(e, error.norm())
    print(f"Norm of Error: {e:<15.6e}")
    return e

def harnessLoop(mockProblems, parameterLists):
    e = -np.inf
    for mp, pl in zip(mockProblems, parameterLists):
        e = max(e, harness(mp, pl))
    return e