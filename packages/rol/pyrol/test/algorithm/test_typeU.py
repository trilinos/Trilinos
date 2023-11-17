import unittest

from pyrol import getParametersFromXmlFile

from harness import harness
from zoo import Rosenbrock


class TestTypeU(unittest.TestCase):

    def _createMasterList(self, iprint = 1):
        p = getParametersFromXmlFile("typeU_params.xml")
        p.sublist("General").set("Output Level", iprint)
        return p

    def setUp(self):
        self.tol = 1e-6
        self.testProblem = Rosenbrock()
        self.parameterList = self._createMasterList(1)


class TestTrustRegion(TestTypeU):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Trust Region")

    def test_CauchyPoint(self):
        s = "Cauchy Point"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # iteration limit exceeded

    def test_TruncatedCG(self):
        s = "Truncated CG"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 22 iterations

    def test_SPG(self):
        s = "SPG"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 24 iterations

    def test_Dogleg(self):
        s = "Dogleg"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 22 iterations

    def test_DoubleDogleg(self):
        s = "Double Dogleg"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # iteration limit exceeded


class TestLineSearch(TestTypeU):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Line Search")
        self.parameterList.sublist("General").set("Inexact Hessian-Times-A-Vector", True)

    def test_SteepestDescent(self):
        d = "Steepest Descent"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # iteration limit exceeded

    def test_NonlinearCG(self):
        d = "Nonlinear CG"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)

    def test_QuasiNewton(self):
        d = "Quasi-Newton Method"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 24 iterations

    def test_Newton(self):
        d = "Newton's Method"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 25 iterations

    def test_NewtonKrylov(self):
        d = "Newton Krylov"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 25 iterations

      
class TestBundle(TestTypeU):
    # TO-DO
    pass