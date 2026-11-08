"""Linear least-squares regression written as a SimOpt problem.

The state u holds the model's predictions and the control z = (a, b) its
parameters. The constraint is u = a*x + b and the objective is
||u - y||^2 / 2. ROL solves the reduced problem in z, and the result is
compared with numpy.linalg.lstsq.
"""
import numpy as np

from pyrol import Problem, Solver, getCout, getParametersFromXmlFile
from pyrol.unsupported import Constraint_SimOpt, Objective_SimOpt, Reduced_Objective_SimOpt
from pyrol.vectors import NumPyVector


class ForwardProblem(Constraint_SimOpt):
    """c(u, z) = u - X z, where the rows of X are (x_i, 1)."""

    def __init__(self, x):
        self.X = np.column_stack([x, np.ones_like(x)])
        super().__init__()

    def value(self, c, u, z, tol):
        c.array[:] = u.array - self.X @ z.array

    def applyJacobian_1(self, jv, v, u, z, tol):
        jv.array[:] = v.array

    def applyJacobian_2(self, jv, v, u, z, tol):
        jv.array[:] = -self.X @ v.array

    def applyInverseJacobian_1(self, ijv, v, u, z, tol):
        ijv.array[:] = v.array

    def applyAdjointJacobian_1(self, ajv, v, u, z, tol):
        ajv.array[:] = v.array

    def applyAdjointJacobian_2(self, ajv, v, u, z, tol):
        ajv.array[:] = -self.X.T @ v.array

    def applyInverseAdjointJacobian_1(self, iajv, v, u, z, tol):
        iajv.array[:] = v.array

    def solve(self, c, u, z, tol):
        u.array[:] = self.X @ z.array
        self.value(c, u, z, tol)


class Loss(Objective_SimOpt):
    """f(u, z) = ||u - y||^2 / 2."""

    def __init__(self, y):
        self.y = y
        super().__init__()

    def value(self, u, z, tol):
        return 0.5 * float(np.sum((u.array - self.y) ** 2))

    def gradient_1(self, g, u, z, tol):
        g.array[:] = u.array - self.y

    def gradient_2(self, g, u, z, tol):
        g.zero()

    def hessVec_11(self, hv, v, u, z, tol):
        hv.array[:] = v.array

    def hessVec_12(self, hv, v, u, z, tol):
        hv.zero()

    def hessVec_21(self, hv, v, u, z, tol):
        hv.zero()

    def hessVec_22(self, hv, v, u, z, tol):
        hv.zero()


def main():
    a, b = 1.0, 0.0  # the "true" parameters from which the data is generated
    n = 100          # the number of data points

    rng = np.random.default_rng(23040901)
    x = rng.uniform(-1, 1, n)
    y = a * x + b + 0.1 * rng.normal(0, 1, n)

    constraint = ForwardProblem(x)
    objective = Loss(y)

    u = NumPyVector(np.zeros(n))   # state
    z = NumPyVector(np.zeros(2))   # control: (a, b)
    l = u.dual().clone()           # adjoint
    g = z.dual().clone()           # gradient with respect to z

    reduced_objective = Reduced_Objective_SimOpt(objective, constraint, u, z, l)
    problem = Problem(reduced_objective, z, g)
    stream = getCout()
    problem.check(True, stream)

    solver = Solver(problem, getParametersFromXmlFile("input.xml"))
    solver.solve(stream)

    target, *_ = np.linalg.lstsq(constraint.X, y, rcond=None)
    np.testing.assert_allclose(z.array, target)
    print('Test passed: Correctly optimized!')


if __name__ == "__main__":
    main()
