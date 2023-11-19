import unittest

from pyrol import getParametersFromXmlFile, getCout

from harness import harness
from zoo import HS32


class TestTypeG(unittest.TestCase):
    
    def setUp(self):
        self.testProblem = HS32()
        self.tol = 1e-8
        # Set up a parameter list.
        p = getParametersFromXmlFile("empty.xml")
        p.sublist("General").set("Output Level", 1)
        self.parameterList = p


class TestAugmentedLagrangian(TestTypeG):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Augmented Lagrangian")

    def test_AugmentedLagrangian(self):
        self.parameterList.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling", False)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 1 iteration
        # test01
        # list.sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Iteration Limit",20);
        # list.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling",false);
        # list.sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);


class TestInteriorPoint(TestTypeG):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Interior Point")

    def test_InteriorPoint(self):
        self.parameterList.sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor", 0.1)
        self.parameterList.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling", False)
        self.parameterList.sublist("Status Test").set("Gradient Tolerance", 1e-9)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < 1e-2)
        # test03
        # list.sublist("Step").sublist("Interior Point").sublist("Subproblem").set("Iteration Limit",200);
        # list.sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor",0.1);
        # list.sublist("Step").sublist("Interior Point").sublist("Subproblem").set("Print History",false);
        # list.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling",false);
        # list.sublist("Status Test").set("Gradient Tolerance",1e-10);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-10);
        # list.sublist("Status Test").set("Step Tolerance",1e-14);
        # list.sublist("Status Test").set("Iteration Limit", 250);
        # list.sublist("Step").set("Type","Trust Region");


class TestMoreauYosida(TestTypeG):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Moreau-Yosida")

    def test_MoreauYosida(self):
        self.parameterList.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling", False)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 5 iterations
        # test02
        # list.sublist("Step").sublist("Moreau-Yosida Penalty").set("Initial Penalty Parameter",1e1);
        # list.sublist("Step").sublist("Moreau-Yosida Penalty").set("Penalty Parameter Growth Factor",1e1);
        # list.sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Iteration Limit",200);
        # list.sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Print History",false);
        # list.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling",false);
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);


class TestStabilizedLCL(TestTypeG):

    def setUp(self):
        super().setUp()
        self.parameterList.sublist("Step").set("Type", "Stabilized LCL")

    def test_StabilizedLCL(self):
        self.parameterList.sublist("Step").sublist("Stabilized LCL").set("Use Default Problem Scaling", False)
        self.parameterList.sublist("Step").sublist("Stabilized LCL").set("Use Default Initial Penalty Parameter", False)
        self.parameterList.sublist("Status Test").set("Gradient Tolerance", 1e-8)
        e = harness(self.testProblem, self.parameterList)
        self.assertTrue(abs(e) < self.tol)  # 8 iterations
        # test04
        # list.sublist("Step").sublist("Stabilized LCL").set("Subproblem Iteration Limit",20);
        # list.sublist("Step").sublist("Stabilized LCL").set("Use Default Problem Scaling",false);
        # list.sublist("Step").sublist("Stabilized LCL").set("Use Default Initial Penalty Parameter",false);
        # list.sublist("Step").sublist("Stabilized LCL").set("Initial Penalty Parameter",10.0);
        # list.sublist("Step").sublist("Stabilized LCL").set("Subproblem Step Type","Trust Region");
        # list.sublist("Status Test").set("Gradient Tolerance",1e-8);
        # list.sublist("Status Test").set("Constraint Tolerance",1e-8);
        # list.sublist("Status Test").set("Step Tolerance",1e-12);
        # list.sublist("Status Test").set("Iteration Limit", 250);