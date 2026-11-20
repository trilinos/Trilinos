import gc
import unittest

import numpy as np

from pyrol import Bounds, Problem, Solver, getCout
from pyrol.vectors import NumPyVector

from helpers import Quadratic, parameters


class TestLifetime(unittest.TestCase):
    """Python objects that only C++ refers to must stay alive."""

    def solve(self, problem):
        solver = Solver(problem, parameters())
        solver.solve(getCout())
        return solver

    def test_objective_only_held_by_problem(self):
        x = NumPyVector(np.zeros(3))
        problem = Problem(Quadratic([1.0, 2.0, 3.0]), x)
        gc.collect()
        self.solve(problem)
        np.testing.assert_allclose(x.array, [1.0, 2.0, 3.0], atol=1e-8)

    def test_problem_built_in_a_function(self):
        def make_problem():
            return Problem(Quadratic([1.0, -1.0]), NumPyVector(np.zeros(2)))

        problem = make_problem()
        gc.collect()
        state = self.solve(problem).getAlgorithmState()
        self.assertLess(state.value, 1e-12)

    def test_bounds_built_from_temporaries(self):
        x = NumPyVector(np.zeros(3))
        problem = Problem(Quadratic([2.0, -2.0, 0.5]), x)
        problem.addBoundConstraint(Bounds(NumPyVector(np.full(3, -1.0)), NumPyVector(np.full(3, 1.0))))
        gc.collect()
        self.solve(problem)
        np.testing.assert_allclose(x.array, [1.0, -1.0, 0.5], atol=1e-6)

    def test_construction_leaves_the_class_alone(self):
        # The add* methods are wrapped once, when PyROL wraps the class.
        before = dict(vars(Problem))
        Problem(Quadratic([1.0]), NumPyVector(np.zeros(1)))
        self.assertEqual(dict(vars(Problem)), before)
        self.assertEqual(Problem.addBoundConstraint.__name__, "addBoundConstraint")

    def test_clones_made_by_cpp(self):
        class TaggedVector(NumPyVector):
            clones = 0

            def clone(self):
                TaggedVector.clones += 1
                return TaggedVector(np.zeros_like(self.array))

        x = TaggedVector(np.full(4, -2.0))
        self.solve(Problem(Quadratic(np.full(4, 3.0)), x))
        gc.collect()
        self.assertGreater(TaggedVector.clones, 0)
        np.testing.assert_allclose(x.array, 3.0, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
