from TorchVectors    import TensorDictVector
from TorchObjectives import LeastSquaresObjective
from Models          import NN

from pyrol import getCout, Problem, Solver
from pyrol.pyrol.Teuchos import ParameterList

import numpy as np

import torch


def generate_data(model):
    n = 10  # number of data points
    torch.manual_seed(24010101)
    x = torch.normal(0, 1, size=(n, model.input_size))
    x.requires_grad = False
    # TO-DO: Consider adding noise to y.
    model.eval()
    with torch.no_grad():
        y = model.forward(x)
    return x, y


def build_parameter_list():
    params = ParameterList()
    params['General'] =  ParameterList()
    params['General']['Output Level'] = 1
    params['Step'] = ParameterList()
    params['Step']['Trust Region'] = ParameterList()
    params['Step']['Trust Region']['Subproblem Solver'] = 'Truncated CG'
    params['Step']['Trust Region']['Initial Radius'] = 1e2
    return params


def check_result(x, alpha):
    target = []
    actual = []
    for i in range(x.dimension()):
        target.append(alpha)
        actual.append(x[i])
    np.testing.assert_allclose(actual, target)
    print('Check passed! Correctly optimized.')


def main():
    torch.set_default_dtype(torch.float64)

    input_size = 2
    model = NN(input_size)
    m = TensorDictVector(model.state_dict())  # model parameters
    alpha = 1.1
    m.setScalar(alpha)
    data = generate_data(model)

    x = m.clone()  # our optimization variable

    objective = LeastSquaresObjective(data, model)
    g = x.clone()

    stream = getCout()

    problem = Problem(objective, x, g)
    problem.checkDerivatives(True, stream)

    params = build_parameter_list()
    solver = Solver(problem, params)
    solver.solve(stream)

    check_result(x, alpha)


if __name__ == '__main__':
    main()
