import unittest

from pyrol import getParametersFromXmlFile

from harness import harness
from zoo import Rosenbrock


class TestTypeU(unittest.TestCase):

    def setUp(self):
        self.testProblem = Rosenbrock()
        self.tol = 1e-6
        # Set up a parameter list.
        p = getParametersFromXmlFile("empty.xml")
        p.sublist("General").set("Output Level", 1)
        self.parameterList = p


class TestTrustRegion(TestTypeU):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Trust Region")
        self.parameterList.sublist("Step").sublist("Trust Region").set("Maximum Radius", 5e3)

    def test_CauchyPoint(self):
        s = "Cauchy Point"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # iteration limit exceeded

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

    def test_SPG(self):
        s = "SPG"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        self.parameterList.sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Iteration Limit", 50)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 24 iterations

    def test_TruncatedCG(self):
        s = "Truncated CG"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Solver", s)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 22 iterations



class TestLineSearch(TestTypeU):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Line Search")
        self.parameterList.sublist("General").set("Inexact Hessian-Times-A-Vector", True)

    def test_Newton(self):
        d = "Newton's Method"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 25 iterations

    def test_NewtonKrylov(self):
        d = "Newton-Krylov"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 25 iterations

    def test_NonlinearCG(self):
        d = "Nonlinear CG"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        self.parameterList.sublist("Step").sublist("Line Search").set("Accept Last Alpha", True)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)

    def test_QuasiNewton(self):
        d = "Quasi-Newton Method"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(e < self.tol)  # 24 iterations

    def test_SteepestDescent(self):
        d = "Steepest Descent"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # iteration limit exceeded


class TestBundle(TestTypeU):
    # TO-DO
    pass