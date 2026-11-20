import unittest

import numpy as np

from pyrol import Problem, Solver, getCout
from pyrol.unsupported import (Constraint_SimOpt, Objective_SimOpt, Reduced_Objective_SimOpt,
                               Vector_SimOpt)
from pyrol.vectors import NumPyVector

from helpers import parameters


class Line(Constraint_SimOpt):
    """c(u, z) = u - X z: the state u holds the predictions of a line z."""

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


class Misfit(Objective_SimOpt):
    """f(u, z) = ||u - y||^2 / 2."""

    def __init__(self, y):
        self.y = np.asarray(y, dtype=float)
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


class TestValueDispatch(unittest.TestCase):
    """value() can be called with (u, z) or with x = Vector_SimOpt(u, z)."""

    def setUp(self):
        self.u = NumPyVector(np.array([1.0, 2.0, 3.0]))
        self.z = NumPyVector(np.array([0.5, 0.0]))

    def test_objective(self):
        objective = Misfit(np.zeros(3))
        self.assertEqual(objective.value(self.u, self.z, 0.0), 7.0)
        self.assertEqual(objective.value(Vector_SimOpt(self.u, self.z), 0.0), 7.0)
        with self.assertRaises(TypeError):
            objective.value(self.u)

    def test_constraint(self):
        constraint = Line(np.array([0.0, 1.0, 2.0]))
        c1, c2 = NumPyVector(np.zeros(3)), NumPyVector(np.zeros(3))
        constraint.value(c1, self.u, self.z, 0.0)
        constraint.value(c2, Vector_SimOpt(self.u, self.z), 0.0)
        np.testing.assert_allclose(c1.array, [1.0, 1.5, 2.0])
        np.testing.assert_allclose(c2.array, c1.array)
        with self.assertRaises(TypeError):
            constraint.value(c1, self.u)


class TestReducedObjective(unittest.TestCase):
    """Least-squares fit of a line, as a reduced SimOpt problem in z."""

    def setUp(self):
        rng = np.random.default_rng(23040901)
        x = rng.uniform(-1, 1, 50)
        self.y = x + 0.1 * rng.normal(0, 1, 50)
        self.constraint = Line(x)
        u = NumPyVector(np.zeros(50))
        self.z = NumPyVector(np.zeros(2))
        self.reduced = Reduced_Objective_SimOpt(Misfit(self.y), self.constraint, u, self.z,
                                                u.dual().clone())

    def test_gradient_check(self):
        d = NumPyVector(np.array([0.3, -0.7]))
        table = self.reduced.checkGradient(self.z, d, False)
        self.assertLess(min(row[3] for row in table), 1e-6)

    def test_solve_matches_least_squares(self):
        problem = Problem(self.reduced, self.z, self.z.dual().clone())
        Solver(problem, parameters()).solve(getCout())
        target, *_ = np.linalg.lstsq(self.constraint.X, self.y, rcond=None)
        np.testing.assert_allclose(self.z.array, target, rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
