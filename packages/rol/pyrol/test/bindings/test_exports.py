import importlib.metadata
import unittest
from pathlib import Path

import pyrol
import pyrol.unsupported


class TestExports(unittest.TestCase):

    def test_supported_classes(self):
        for name in ["Bounds", "Constraint", "LinearConstraint", "LinearOperator",
                     "Objective", "Problem", "Solver", "Vector"]:
            with self.subTest(name=name):
                cls = getattr(pyrol, name)
                forScalar = getattr(pyrol, name + "_forScalar")
                self.assertIs(forScalar(), cls)
                self.assertIs(forScalar("double"), cls)
                self.assertEqual(cls.__name__, name + "_double_t")
                with self.assertRaises(ValueError):
                    forScalar("float")

    def test_supported_functions(self):
        for name in ["getCout", "getParametersFromXmlFile", "getParametersFromYamlFile"]:
            with self.subTest(name=name):
                self.assertTrue(callable(getattr(pyrol, name)))

    def test_unsupported_classes(self):
        for name in ["Constraint_SimOpt", "Objective_SimOpt", "PartitionedVector",
                     "Reduced_Objective_SimOpt", "Vector_SimOpt"]:
            with self.subTest(name=name):
                self.assertTrue(hasattr(pyrol.unsupported, name))
                self.assertFalse(hasattr(pyrol, name))

    def test_vectors_star_import(self):
        namespace = {}
        exec("from pyrol.vectors import *", namespace)
        self.assertIn("NumPyVector", namespace)

    def test_version(self):
        self.assertEqual(pyrol.version(), "PyROL version: " + pyrol.__version__)
        # Installed with pip, the version is the package's. Built with
        # Trilinos, there is no package metadata next to pyrol, and the
        # version is Trilinos'.
        site = str(Path(pyrol.__file__).resolve().parents[1])
        for dist in importlib.metadata.distributions(path=[site]):
            if dist.metadata["Name"] == "rol-python":
                self.assertEqual(pyrol.__version__, dist.version)
                break
        else:
            self.assertRegex(pyrol.__version__, r"^\d+\.\d+")


if __name__ == "__main__":
    unittest.main()
