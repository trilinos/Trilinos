import unittest

import numpy as np

from pyrol import LinearConstraint, LinearOperator, Problem, Solver, getCout
from pyrol.vectors import NumPyVector

from helpers import Quadratic, parameters


class RowSum(LinearOperator):
    """The 1 x n matrix of ones."""

    def apply(self, Hv, v, tol):
        Hv.array[0] = np.sum(v.array)

    def applyAdjoint(self, Hv, v, tol):
        Hv.array[:] = v.array[0]


class TestProblem(unittest.TestCase):

    def test_linear_equality_constraint(self):
        # Minimize ||x - c||^2 / 2 subject to sum(x) = 1, with the constraint
        # A x + b = 0 built from a Python linear operator A.
        center = np.array([1.0, 2.0, 3.0])
        x = NumPyVector(np.zeros(3))
        problem = Problem(Quadratic(center), x)
        constraint = LinearConstraint(RowSum(), NumPyVector(np.array([-1.0])))
        problem.addLinearConstraint("sum", constraint, NumPyVector(np.zeros(1)))
        Solver(problem, parameters("Augmented Lagrangian")).solve(getCout())
        np.testing.assert_allclose(x.array, center - (center.sum() - 1) / 3, atol=1e-6)

    def test_solver_state(self):
        x = NumPyVector(np.full(2, 5.0))
        solver = Solver(Problem(Quadratic([1.0, 1.0]), x), parameters())
        solver.solve(getCout())
        state = solver.getAlgorithmState()
        self.assertGreater(state.iter, 0)
        self.assertGreater(state.nfval, 0)
        self.assertLess(state.value, 1e-12)
        self.assertLess(state.gnorm, 1e-6)


if __name__ == "__main__":
    unittest.main()
