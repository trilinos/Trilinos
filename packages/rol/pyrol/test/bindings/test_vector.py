import unittest

import numpy as np

from pyrol import Vector
from pyrol.unsupported import (Constraint_SimOpt, Objective_SimOpt, PartitionedVector,
                               Reduced_Objective_SimOpt, Vector_SimOpt)
from pyrol.vectors import NumPyVector


class CountingVector(Vector):
    """Minimal vector that relies on the default ROL::Vector::axpy."""

    clones = 0

    def __init__(self, values):
        self.values = np.array(values, dtype=float)
        super().__init__()

    def plus(self, x):
        self.values += x.values

    def scale(self, alpha):
        self.values *= alpha

    def dot(self, x):
        return float(np.dot(self.values, x.values))

    def norm(self):
        return float(np.sqrt(self.dot(self)))

    def clone(self):
        CountingVector.clones += 1
        return CountingVector(np.zeros_like(self.values))


class DualVector(NumPyVector):
    """Riesz map x -> 2x; dual() returns a new vector on every call."""

    def dual(self):
        return type(self)(2.0 * self.array)


class BasisVector(NumPyVector):
    """basis() returns a new vector on every call."""

    def basis(self, i):
        b = type(self)(np.zeros_like(self.array))
        b.array[i] = 1.0
        return b


class TestVectorClone(unittest.TestCase):

    def test_cpp_clone_calls_python_clone_once(self):
        y = CountingVector([1.0, 2.0])
        x = CountingVector([3.0, 4.0])
        CountingVector.clones = 0
        y.axpy(2.0, x)  # ROL::Vector::axpy clones x exactly once
        self.assertEqual(CountingVector.clones, 1)
        np.testing.assert_allclose(y.values, [7.0, 10.0])


class TestOverridesCalledFromCpp(unittest.TestCase):
    """C++ keeps using what dual() and basis() return after the Python call."""

    def test_dual_returning_new_vector(self):
        # The constructor calls dual().clone() on its vectors; this used to
        # crash because C++ referred to a vector that Python had freed.
        u, z = DualVector(np.ones(3)), DualVector(np.ones(2))
        Reduced_Objective_SimOpt(Objective_SimOpt(), Constraint_SimOpt(), u, z, u.dual().clone())

    def test_dual_values_seen_by_cpp(self):
        # Vector::apply computes x.dot(y.dual()) in C++, and Vector_SimOpt::dual
        # assembles the duals of its parts.
        x = Vector_SimOpt(DualVector(np.array([1.0, 2.0])), DualVector(np.array([3.0])))
        y = Vector_SimOpt(DualVector(np.array([4.0, 5.0])), DualVector(np.array([6.0])))
        self.assertAlmostEqual(x.apply(y), 2 * (1 * 4 + 2 * 5 + 3 * 6))

    def test_basis_returning_new_vector(self):
        # PartitionedVector::basis copies the basis vectors of its parts.
        x = PartitionedVector.create(BasisVector(np.zeros(2)), 2)
        e = x.basis(3)
        np.testing.assert_allclose(e.get(0).array, [0.0, 0.0])
        np.testing.assert_allclose(e.get(1).array, [0.0, 1.0])


class TestStdVector(unittest.TestCase):
    """std::vector results, e.g. of checkVector, behave like Python sequences."""

    def test_sequence_protocol(self):
        rng = np.random.default_rng(0)
        x, y, z = (NumPyVector(rng.standard_normal(3)) for _ in range(3))
        errors = x.checkVector(y, z, False)
        self.assertGreater(len(errors), 0)
        self.assertEqual(len(list(errors)), len(errors))
        self.assertEqual(errors[-1], errors[len(errors) - 1])
        self.assertEqual(errors[1:3], [errors[1], errors[2]])
        with self.assertRaises(IndexError):
            errors[len(errors)]


if __name__ == "__main__":
    unittest.main()
