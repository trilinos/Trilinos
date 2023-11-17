import unittest

from pyrol import getParametersFromXmlFile

# TO-DO: clean up the algorithm names
from pyrol import getCout

from harness import harness
from zoo import HS42

class TestTypeE(unittest.TestCase):

    def _createMasterList(self, iprint = 1):
        p = getParametersFromXmlFile("typeE_parameters.xml")
        p.sublist("General").set("Output Level", iprint)
        return p

    def setUp(self):
        self.tol = 1e-8
        self.mockProblem = HS42()
        self.parameterList = self._createMasterList(1)

    def test_AugmentedLagrangian(self):
        self.parameterList.sublist("Step").set("Type", "Augmented Lagrangian")
        e = harness(self.mockProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)

    def test_CompositeStep(self):
        self.parameterList.sublist("Step").set("Type", "Composite Step")
        e = harness(self.mockProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)
    
    def test_Fletcher(self):
        self.parameterList.sublist("Step").set("Type", "Fletcher")
        e = harness(self.mockProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)

    def test_StabilizedLCL(self):
        self.parameterList.sublist("Step").set("Type", "Stabilized LCL")
        e = harness(self.mockProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)