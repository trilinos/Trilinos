import numpy as np

from pyrol import Objective, Problem, Solver, getCout, getParametersFromXmlFile
from pyrol.vectors import NumPyVector

class MyObjective(Objective):
    def value(self, x, tol):
        return 0

def makeProblem():
    objective = MyObjective()
    x = NumPyVector(np.zeros(1))
    p = Problem(objective, x)
    return p
    
def main():
    p = makeProblem()
    stream = getCout()
    list = getParametersFromXmlFile("input.xml")
    solver = Solver(p, list)
    solver.solve(stream)

if __name__ == "__main__":
    main()