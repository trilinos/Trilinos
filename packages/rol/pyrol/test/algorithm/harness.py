import numpy as np

from pyrol import getCout, Solver
from pyrol.vectors import NumPyVector

def harness(testProblem, parameterList):
    stream = getCout()
    problem, x, solutions = testProblem.get()

    # TO-DO: print problem name

    solver = Solver(problem, parameterList)
    solver.solve(stream)

    error = x.clone()
    e = np.inf
    for s in solutions:
        error[:] = x[:] - s[:]
        e = min(e, error.norm())
    print(f"Norm of Error: {e:<15.6e}")
    return e

def harnessLoop(testProblems, parameterLists):
    e = -np.inf
    for mp, pl in zip(testProblems, parameterLists):
        e = max(e, harness(mp, pl))
    return e