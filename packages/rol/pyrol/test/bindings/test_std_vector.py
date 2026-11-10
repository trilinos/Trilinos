import gc
import unittest

import numpy as np

from pyrol import Objective, Problem, Solver, getCout
from pyrol.pyrol import std
from pyrol.unsupported import PartitionedVector, TimeStamp
from pyrol.vectors import NumPyVector

from helpers import parameters


class TestNumberVectors(unittest.TestCase):
    """Python sequences and NumPy arrays where ROL takes a std::vector<double>."""

    def test_converted_from_python(self):
        ts = TimeStamp()
        for value in ([0.0, 0.5], (0.0, 0.5), np.array([0.0, 0.5])):
            with self.subTest(value=value):
                ts.t = value
                self.assertEqual(list(ts.t), [0.0, 0.5])

    def test_rejects_non_numbers(self):
        ts = TimeStamp()
        with self.assertRaises(TypeError):
            ts.t = ["a", "b"]

    def test_sequence(self):
        v = std.vector_double_t([1.0, 2.0, 3.0])
        self.assertEqual(len(v), 3)
        self.assertEqual(list(v), [1.0, 2.0, 3.0])
        self.assertEqual(v[-1], 3.0)
        self.assertEqual(list(v[1:]), [2.0, 3.0])
        with self.assertRaises(IndexError):
            v[3]
        v.append(4.0)
        self.assertEqual(v, [1.0, 2.0, 3.0, 4.0])

    def test_cpp_style_methods(self):
        v = std.vector_double_t(2, 1.5)
        v.push_back(3.0)
        self.assertEqual((v.size(), v.at(2), v.front(), v.back()), (3, 3.0, 1.5, 3.0))
        v.pop_back()
        v.resize(1)
        self.assertEqual(list(v), [1.5])
        with self.assertRaises(IndexError):
            v.at(1)

    def test_numpy_view(self):
        ts = TimeStamp()
        ts.t = [0.0, 0.5]
        view = np.asarray(ts.t)
        view[1] = 1.0
        self.assertEqual(ts.t[1], 1.0)


class TestVectorsOfObjects(unittest.TestCase):
    """Python lists where ROL takes a std::vector of ROL objects."""

    def test_partitioned_vector_from_list(self):
        x = PartitionedVector([NumPyVector(np.array([1.0, 2.0])), NumPyVector(np.array([3.0]))])
        self.assertEqual((x.numVectors(), x.dimension()), (2, 3))
        self.assertAlmostEqual(x.dot(x), 14.0)

    def test_list_elements_kept_alive(self):
        x = PartitionedVector([NumPyVector(np.array([1.0, 2.0])), NumPyVector(np.array([3.0]))])
        gc.collect()
        self.assertIsInstance(x.get(1), NumPyVector)
        y = x.clone()
        y.set(x)
        self.assertAlmostEqual(y.dot(x), 14.0)

    def test_solve_over_partitioned_vector(self):
        targets = [np.array([1.0, -2.0]), np.array([3.0])]

        class Distance(Objective):
            """f(x) = sum_i ||x_i - target_i||^2 / 2 over the parts of x."""

            def value(self, x, tol):
                return 0.5 * sum(float(np.sum((x.get(i).array - t) ** 2)) for i, t in enumerate(targets))

            def gradient(self, g, x, tol):
                for i, t in enumerate(targets):
                    g.get(i).array[:] = x.get(i).array - t

        x = PartitionedVector([NumPyVector(np.zeros(2)), NumPyVector(np.zeros(1))])
        Solver(Problem(Distance(), x), parameters()).solve(getCout())
        for i, t in enumerate(targets):
            np.testing.assert_allclose(x.get(i).array, t, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
