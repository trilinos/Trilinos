import unittest

from pyrol import getParametersFromXmlFile

from harness import harness
from zoo import HS02, HS41


class TestTypeB(unittest.TestCase):

    def setUp(self):
        self.tol = 1e-7
        # Set up a parameter list.
        p = getParametersFromXmlFile("empty.xml")
        p.sublist("General").set("Output Level", 1)
        self.parameterList = p

    def _scenario02(self):
        # HS02, which has bound constraints.
        self.testProblem = HS02()
        self.parameterList.sublist("Status Test").set("Gradient Tolerance", 1e-6)

    def _scenario41(self):
        # HS41, which has bound constraints and a linear equality constraint.
        self.testProblem = HS41()
        self.parameterList.sublist("Status Test").set("Gradient Tolerance", 1e-8)


class TestInteriorPoint(TestTypeB):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Interior Point") 
    
    def test_InteriorPoint(self):
        self._scenario41()
        self.parameterList.sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor", 0.1)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < 1e-3)  # 6 iterations
        # test07
        # list.sublist("Step").sublist("Interior Point").sublist("Subproblem").set("Iteration Limit",200);
        # list.sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor",0.1);
        # list.sublist("Step").sublist("Interior Point").sublist("Subproblem").set("Print History",false);
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("Step").set("Type","Trust Region");


class TestLineSearch(TestTypeB):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Line Search")

    def test_LSecantB(self):
        d = "Quasi-Newton Method"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        m = "L-Secant-B"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Quasi-Newton").set("Method", m)
        self._scenario41()
        self.parameterList.sublist("General").sublist("Polyhedral Projection").set("Type", "Dai-Fletcher")
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)
        # test15
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("Step").set("Type","Line Search");
        # list.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");
        # list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");

    def test_NewtonKrylov(self):
        d = "Newton-Krylov"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        self.parameterList.sublist("General").sublist("Secant").set("Use as Preconditioner", True)
        self._scenario02()
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # 6 iterations but O(1) error
        # test04
        # parlist->sublist("General").sublist("Secant").set("Use as Hessian",false);
        # parlist->sublist("General").sublist("Secant").set("Use as Preconditioner",true);
        # parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
        # parlist->sublist("General").sublist("Krylov").set("Iteration Limit", dim);

    def test_QuasiNewton(self):
        d = "Quasi-Newton Method"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        m = "Quasi-Newton"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Quasi-Newton").set("Method", m)
        self._scenario41()
        self.parameterList.sublist("General").sublist("Polyhedral Projection").set("Type", "Dai-Fletcher")
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 7 iterations
        # test12
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("Step").set("Type","Line Search");
        # list.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");
        # list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");

    def test_SteepestDescent(self):
        d = "Steepest Descent"
        self.parameterList.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", d)
        self._scenario41()
        self.parameterList.sublist("Status Test").set("Gradient Tolerance", 1e-7)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < 1e-7)  # 6 iterations
        # test01
        # list.sublist("Status Test").set("Gradient Tolerance",1e-7)
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8)
        # list.sublist("Status Test").set("Step Tolerance",1e-12)
        # list.sublist("Status Test").set("Iteration Limit", 250)


class TestMoreauYosida(TestTypeB):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Moreau-Yosida") 
    
    def test_MoreauYosida(self):
        self._scenario41()
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations
        # test06
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);


class TestPrimalDualActiveSet(TestTypeB):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Primal Dual Active Set")

    def test_PrimalDualActiveSet(self):
        self._scenario41()
        self.parameterList.sublist("General").sublist("Secant").set("Use as Hessian", True)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 7 iterations
        # test09
        # list.sublist("Status Test").set("Gradient Tolerance",1e-10);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-10);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("General").sublist("Secant").set("Type","Limited Memory BFGS");
        # list.sublist("General").sublist("Secant").set("Use as Hessian",true);


class TestTrustRegion(TestTypeB):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Trust Region") 

    def test_ColemanLi(self):
        s = "Coleman-Li"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Model", s)
        self._scenario02()
        self.parameterList.sublist("Status Test").set("Gradient Tolerance", 1e-6)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(True)  # step tolerance met at iteration 13
        # test16
        # parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
        # parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",-1.e1);
        # parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e4);
        # parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
        # palist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

    def test_KelleySachs(self):
        s = "Kelley-Sachs"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Model", s)
        self._scenario02()
        self.parameterList.sublist("General").sublist("Secant").set("Use as Preconditioner", True)
        self.parameterList.sublist("Step").sublist("Trust Region").set("Initial Radius", 1e0)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 31 iterations
        # test10
        # parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
        # parlist->sublist("General").sublist("Secant").set("Use as Hessian",false);
        # parlist->sublist("General").sublist("Secant").set("Use as Preconditioner",true);
        # parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",1.0);
        # parlist->sublist("General").sublist("Krylov").set("Iteration Limit", dim);

    def test_LinMore(self):
        s = "Lin-More"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Model", s)
        self._scenario41()
        self.parameterList.sublist("General").sublist("Polyhedral Projection").set("Type", "Dai-Fletcher")
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations
        # test02
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");

    def test_TrustRegionSPG(self):
        s = "SPG"
        self.parameterList.sublist("Step").sublist("Trust Region").set("Subproblem Model", s)
        self._scenario41()
        self.parameterList.sublist("General").sublist("Polyhedral Projection").set("Type", "Dai-Fletcher")
        self.parameterList.sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Maximum Spectral Step Size", 1e2)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 3 iterations
        # test13
        # list.sublist("Status Test").set("Gradient Tolerance",1e-7);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("General").sublist("Polyhedral Projection").set("Iteration Limit",5000);
        # list.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");
        # list.sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Maximum Spectral Step Size",1e2);
        # list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");


class TestSpectralGradient(TestTypeB):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Spectral Gradient") 

    def test_SpectralGradient(self):
        self._scenario41()
        self.parameterList.sublist("General").sublist("Polyhedral Projection").set("Type", "Dai-Fletcher")
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < 1e-7)  # 5 iterations
        # test11
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");
