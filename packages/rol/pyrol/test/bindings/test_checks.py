import unittest

import numpy as np

from pyrol import Constraint, Objective, Problem
from pyrol.vectors import NumPyVector


class Quartic(Objective):
    """f(x) = sum(x^4) / 4 + ||x||^2 / 2."""

    def value(self, x, tol):
        return float(np.sum(x.array ** 4) / 4 + np.sum(x.array ** 2) / 2)

    def gradient(self, g, x, tol):
        g.array[:] = x.array ** 3 + x.array

    def hessVec(self, hv, v, x, tol):
        hv.array[:] = (3 * x.array ** 2 + 1) * v.array


class Sphere(Constraint):
    """c(x) = ||x||^2 - 1."""

    def value(self, c, x, tol):
        c.array[0] = x.array @ x.array - 1

    def applyJacobian(self, jv, v, x, tol):
        jv.array[0] = 2 * x.array @ v.array

    def applyAdjointJacobian(self, ajv, v, x, tol):
        ajv.array[:] = 2 * v.array[0] * x.array


def smallest_error(table):
    """Smallest finite-difference error in a table returned by ROL's checks."""
    return min(row[3] for row in table)


class TestDerivativeChecks(unittest.TestCase):
    """ROL's finite-difference checks of Python-defined functions."""

    def setUp(self):
        rng = np.random.default_rng(1)
        self.x = NumPyVector(rng.standard_normal(4))
        self.d = NumPyVector(rng.standard_normal(4))

    def test_gradient(self):
        self.assertLess(smallest_error(Quartic().checkGradient(self.x, self.d, False)), 1e-6)

    def test_hessian_vector_product(self):
        self.assertLess(smallest_error(Quartic().checkHessVec(self.x, self.d, False)), 1e-6)

    def test_wrong_gradient_is_detected(self):
        class WrongGradient(Quartic):
            def gradient(self, g, x, tol):
                g.array[:] = x.array  # missing the x^3 term

        self.assertGreater(smallest_error(WrongGradient().checkGradient(self.x, self.d, False)), 1e-3)

    def test_constraint_jacobian(self):
        jv = NumPyVector(np.zeros(1))
        self.assertLess(smallest_error(Sphere().checkApplyJacobian(self.x, self.d, jv, False)), 1e-6)

    def test_constraint_adjoint_consistency(self):
        w = NumPyVector(np.array([0.7]))
        self.assertLess(Sphere().checkAdjointConsistencyJacobian(w, self.d, self.x), 1e-10)

    def test_problem_check(self):
        problem = Problem(Quartic(), self.x)
        problem.addConstraint("sphere", Sphere(), NumPyVector(np.zeros(1)))
        problem.check(False)


if __name__ == "__main__":
    unittest.main()
