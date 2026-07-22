from TorchVectors    import TensorVector
from TorchObjectives import TorchObjective

from pyrol import getCout, Objective, Problem, Solver
from pyrol.vectors import NumPyVector

from pyrol.pyrol.Teuchos import ParameterList

import numpy as np
import time
import torch


class RosenbrockObjective(TorchObjective):

    def __init__(self):
        super().__init__()
        self.alpha = 100

    def torch_value(self, x):
        # return torch.sum(self.alpha*(x[::2]**2 - x[1::2])**2 + (x[::2] - 1)**2)
        return torch.sum(self.alpha*(x[:-1]**2 - x[1:])**2 + (x[:-1] - 1)**2)


def build_parameter_list():
    params = ParameterList()
    params['General'] =  ParameterList()
    params['General']['Output Level'] = 1
    params['Step'] = ParameterList()
    params['Step']['Trust Region'] = ParameterList()
    params['Step']['Trust Region']['Subproblem Solver'] = 'Truncated CG'
    params['Status Test'] = ParameterList()
    params['Status Test']['Iteration Limit'] = 10000

    return params


def main():

    torch.set_default_dtype(torch.float64)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    #device = torch.device('cpu')

    start = time.time()
    n = int(1e2)
    print(device)
    x = torch.empty(n, requires_grad=False, device=device)
    x[ ::2] = -1.2
    x[1::2] =  1.0
    x = TensorVector(x)

    objective = RosenbrockObjective()
    g = x.clone()

    stream = getCout()

    problem = Problem(objective, x, g)
    # problem.checkDerivatives(True, stream)

    params = build_parameter_list()
    solver = Solver(problem, params)
    solver.solve(stream)
    print(f"Solve time: {time.time() - start}\n")

    print(g.torch_object.device)

if __name__ == "__main__":
    main()
