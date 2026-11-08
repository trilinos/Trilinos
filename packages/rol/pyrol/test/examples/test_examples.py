"""Run the examples, so that they keep working."""
import importlib.util
import os
import runpy
import sys
import unittest
from pathlib import Path

import numpy as np

EXAMPLES = Path(__file__).resolve().parents[2] / "example"


def has_module(name):
    return importlib.util.find_spec(name) is not None


class TestExamples(unittest.TestCase):

    def run_example(self, relative_path, run_name="__main__"):
        """Run an example from its own directory, as a user would."""
        path = EXAMPLES / relative_path
        cwd, sys_path = os.getcwd(), list(sys.path)
        os.chdir(path.parent)
        sys.path.insert(0, str(path.parent))
        try:
            return runpy.run_path(str(path), run_name=run_name)
        finally:
            os.chdir(cwd)
            sys.path[:] = sys_path

    def test_rosenbrock(self):
        self.run_example("rosenbrock.py")

    def test_simopt_regression(self):
        self.run_example("simopt/regression.py")

    def test_dynamic_interface(self):
        example = self.run_example("dynamic_interface/forward_euler_1d.py", run_name="example")
        constraint, u0, z, timestamps = example["optimize"](100)
        self.assertLess(np.linalg.norm(example["as_numpy"](z), np.inf), 1e-2)

    @unittest.skipUnless(has_module("torch"), "requires PyTorch")
    def test_pytorch(self):
        for name in ["rosenbrock_torch.py", "tensor_example.py", "nn_example.py"]:
            with self.subTest(example=name):
                self.run_example(f"pytorch/{name}")


if __name__ == "__main__":
    unittest.main()
