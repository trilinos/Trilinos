import unittest

from pyrol import getParametersFromXmlFile

from harness import harness
from zoo import HS42


class TestTypeE(unittest.TestCase):

    def setUp(self):
        self.testProblem = HS42()
        self.tol = 1e-8
        # Set up a parameter list.
        p = getParametersFromXmlFile("empty.xml")
        p.sublist("General").set("Output Level", 1)
        self.parameterList = p
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", "Truncated CG")


class TestAugmentedLagrangian(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Augmented Lagrangian")

    def test_AugmentedLagrangian(self):
        self.parameterList.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Initial Penalty Parameter", False)
        self.parameterList.sublist("Step").sublist("Augmented Lagrangian").set("Initial Penalty Parameter", 5e2)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations
        # test01


class TestCompositeStep(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Composite Step")

    def test_CompositeStep(self):
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations
        # test03


class TestFletcher(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Fletcher")

    def test_Fletcher(self):
        self.parameterList.sublist("Step").sublist("Fletcher").set("Penalty Parameter", 1e2)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 1 iteration
        # test02


class TestStabilizedLCL(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Stabilized LCL")

    def test_StabilizedLCL(self):
        self.parameterList.sublist("Step").sublist("Stabilized LCL").set("Use Default Initial Penalty Parameter", False)
        self.parameterList.sublist("Step").sublist("Stabilized LCL").set("Use Default Problem Scaling", False)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 5 iterations
        # test04