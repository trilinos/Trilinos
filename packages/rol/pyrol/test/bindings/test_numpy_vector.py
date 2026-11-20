import unittest

import numpy as np

from pyrol.pyrol import ROL
from pyrol.unsupported import PartitionedVector
from pyrol.vectors import NumPyVector


class TestNumPyVector(unittest.TestCase):

    def setUp(self):
        self.rng = np.random.default_rng(0)

    def random(self, n=6):
        return NumPyVector(self.rng.standard_normal(n))

    def test_rol_vector_checks(self):
        # ROL's own consistency checks of the vector-space operations.
        x, y, z = self.random(), self.random(), self.random()
        errors = list(x.checkVector(y, z, False))
        self.assertGreater(len(errors), 0)
        self.assertLess(max(abs(e) for e in errors), 1e-10)

    def test_arithmetic(self):
        a, b = self.rng.standard_normal(5), self.rng.standard_normal(5)
        x, y = NumPyVector(a.copy()), NumPyVector(b.copy())
        self.assertEqual(x.dimension(), 5)
        self.assertAlmostEqual(x.dot(y), a @ b)
        self.assertAlmostEqual(x.norm(), np.linalg.norm(a))
        x.axpy(2.0, y)
        np.testing.assert_allclose(x.array, a + 2 * b)
        x.scale(0.5)
        np.testing.assert_allclose(x.array, 0.5 * a + b)
        x.plus(y)
        np.testing.assert_allclose(x.array, 0.5 * a + 2 * b)
        x.setScalar(3.0)
        np.testing.assert_allclose(x.array, 3.0)
        x.zero()
        np.testing.assert_allclose(x.array, 0.0)

    def test_indexing(self):
        x = NumPyVector(np.arange(4.0))
        self.assertEqual(x[1], 1.0)
        x[0:2] = [5.0, 6.0]
        np.testing.assert_allclose(x[:], [5.0, 6.0, 2.0, 3.0])

    def test_clone_and_basis(self):
        x = NumPyVector(np.arange(3.0))
        c = x.clone()
        self.assertIsInstance(c, NumPyVector)
        np.testing.assert_allclose(c.array, 0.0)
        np.testing.assert_allclose(x.basis(1).array, [0.0, 1.0, 0.0])

    def test_elementwise_operations(self):
        x = NumPyVector(np.array([-1.0, 2.0, -3.0]))
        y = NumPyVector(np.array([2.0, 2.0, 2.0]))
        self.assertEqual(x.reduce(ROL.Elementwise.ReductionMin_double_t()), -3.0)
        self.assertEqual(x.reduce(ROL.Elementwise.ReductionMax_double_t()), 2.0)
        self.assertEqual(x.reduce(ROL.Elementwise.ReductionSum_double_t()), -2.0)
        x.applyUnary(ROL.Elementwise.AbsoluteValue_double_t())
        np.testing.assert_allclose(x.array, [1.0, 2.0, 3.0])
        x.applyBinary(ROL.Elementwise.Multiply_double_t(), y)
        np.testing.assert_allclose(x.array, [2.0, 4.0, 6.0])

    def test_operations_called_from_cpp(self):
        # PartitionedVector (C++) forwards these to its NumPyVector parts.
        x = PartitionedVector.create(NumPyVector(np.zeros(2)), 2)
        x.get(0).array[:] = [1.0, -2.0]
        x.get(1).array[:] = [3.0, -4.0]
        self.assertEqual(x.dimension(), 4)
        self.assertAlmostEqual(x.norm(), np.sqrt(30.0))
        self.assertEqual(x.reduce(ROL.Elementwise.ReductionMax_double_t()), 3.0)
        x.applyUnary(ROL.Elementwise.AbsoluteValue_double_t())
        np.testing.assert_allclose(x.get(1).array, [3.0, 4.0])
        y = x.clone()
        y.set(x)
        y.scale(2.0)
        np.testing.assert_allclose(y.get(0).array, [2.0, 4.0])
        np.testing.assert_allclose(x.get(0).array, [1.0, 2.0])


if __name__ == "__main__":
    unittest.main()
