import unittest

from pyrol import getParametersFromXmlFile

from harness import harness
from zoo import HS42


class TestTypeE(unittest.TestCase):

    def _createMasterList(self, iprint = 1):
        p = getParametersFromXmlFile("typeE_params.xml")
        p.sublist("General").set("Output Level", iprint)
        return p

    def setUp(self):
        self.tol = 1e-8
        self.testProblem = HS42()
        self.parameterList = self._createMasterList(1)


class TestAugmentedLagrangian(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Augmented Lagrangian")

    def test_AugmentedLagrangian(self):
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations


class TestCompositeStep(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Composite Step")

    def test_CompositeStep(self):
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations


class TestFletcher(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Fletcher")

    def test_Fletcher(self):
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 1 iteration


class TestStabilizedLCL(TestTypeE):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Stabilized LCL")

    def test_StabilizedLCL(self):
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 5 iterations