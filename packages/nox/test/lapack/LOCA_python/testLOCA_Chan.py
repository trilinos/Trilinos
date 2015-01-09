#! /usr/bin/env python

import setpath
import unittest

import NOX
import LOCA

def runTests(verbose):

    # Our TestCase for the Chan Problem
    # It is declared here so it can access the symbols below
    class ChanTestCase(unittest.TestCase):
        "TestCase for the Chan Problem"

        def testContinuation(self):
            "Test the continuation run finishes correctly"

            self.assertEqual(status, LOCA.Iterator.Finished)

        def testNumSteps(self):
            "Test the number of continuation steps"
            numSteps = stepper.getStepNumber()
            self.assertEqual(numSteps, 32)

        def testNumFailedSteps(self):
            "Test the number of failed continuation steps"
            numFailedSteps = stepper.getNumFailedSteps()
            self.assertEqual(numFailedSteps, 0)

        def testFinalValue(self):
            "Test the final value of the continuation parameter"
            p = grp.getParams()
            alpha_final = p.getValue("alpha")
            self.assertAlmostEqual(alpha_final, 5.0, 10)

        def testNorm(self):
            "Test the norm of the final solution"
            norm_x = finalSolution.norm()
            self.assertAlmostEqual(norm_x, 203.1991024, 6)

    # Set up constants
    n = 100
    alpha = 0.0
    beta = 0.0
    scale = 1.0
    maxNewtonIters = 20

    alpha = alpha / scale

    # Create problem interface
    chan = LOCA.Chan.ProblemInterface(n, alpha, beta, scale)

    # Set up parameter vector
    p = LOCA.ParameterVector()
    p.addParameter("alpha",alpha)
    p.addParameter("beta",beta)
    p.addParameter("scale",scale)

    # Create group
    grp = LOCA.LAPACK.Group(chan)
    grp.setParams(p)

    # Set up parameter lists
    paramList = NOX.Parameter.List()

    locaParamsList = paramList.sublist("LOCA")

    # Create the stepper sublist and set the stepper parameters
    stepperList = locaParamsList.sublist("Stepper")
    stepperList.setParameter("Continuation Method", "Arc Length")
    stepperList.setParameter("Continuation Parameter", "alpha")
    stepperList.setParameter("Initial Value", alpha)
    stepperList.setParameter("Max Value", 5.0/scale)
    stepperList.setParameter("Min Value", 0.0/scale)
    stepperList.setParameter("Max Steps", 50)
    stepperList.setParameter("Max Nonlinear Iterations", maxNewtonIters)
    stepperList.setParameter("Enable Arc Length Scaling", True)
    stepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5)
    stepperList.setParameter("Max Arc Length Parameter Contribution", 0.7)
    stepperList.setParameter("Initial Scale Factor", 1.0)
    stepperList.setParameter("Min Scale Factor", 1.0e-8)
    stepperList.setParameter("Enable Tangent Factor Step Size Scaling",True)
    stepperList.setParameter("Min Tangent Factor", -1.0)
    stepperList.setParameter("Tangent Factor Exponent",1.0)
    stepperList.setParameter("Compute Eigenvalues",False)

    # Create bifurcation sublist
    bifurcationList = locaParamsList.sublist("Bifurcation")
    bifurcationList.setParameter("Method", "None")

    # Create predictor sublist
    predictorList = locaParamsList.sublist("Predictor")
    predictorList.setParameter("Method", "Tangent")
    firstStepPredictorList = predictorList.sublist("First Step Predictor");
    firstStepPredictorList.setParameter("Method", "Constant")

    # Create step size sublist
    stepSizeList = locaParamsList.sublist("Step Size")
    stepSizeList.setParameter("Method", "Adaptive")
    stepSizeList.setParameter("Initial Step Size", 0.1/scale)
    stepSizeList.setParameter("Min Step Size", 1.0e-3/scale)
    stepSizeList.setParameter("Max Step Size", 10.0/scale)
    stepSizeList.setParameter("Aggressiveness", 0.5)
    stepSizeList.setParameter("Failed Step Reduction Factor", 0.5)
    stepSizeList.setParameter("Successful Step Increase Factor", 1.26)

    # Set the LOCA Utilities
    locaUtilsList = locaParamsList.sublist("Utilities")
    if (verbose):
        locaUtilsList.setParameter("Output Information",
                                   LOCA.Utils.Error +
                                   LOCA.Utils.Warning +
                                   LOCA.Utils.StepperIteration +
                                   LOCA.Utils.StepperDetails +
                                   LOCA.Utils.Solver +
                                   LOCA.Utils.Parameters +
                                   LOCA.Utils.SolverDetails)
    else:
        locaUtilsList.setParameter("Output Information",LOCA.Utils.Error)

    # Create the "Solver" parameters sublist to be used with NOX Solvers
    nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based")

    nlPrintParams = nlParams.sublist("Printing");
    if (verbose):
        nlPrintParams.setParameter("Output Information",
                                   NOX.Utils.Error +
                                   NOX.Utils.Details +
                                   NOX.Utils.OuterIteration +
                                   NOX.Utils.InnerIteration +
                                   NOX.Utils.Warning +
                                   NOX.Utils.TestDetails)
    else:
         nlPrintParams.setParameter("Output Information", NOX.Utils.Error)

    # Set up the status tests
    normF = NOX.StatusTest.NormF(1.0e-8)
    maxIters = NOX.StatusTest.MaxIters(maxNewtonIters)
    comboOR = NOX.StatusTest.Combo(NOX.StatusTest.Combo.OR,
                                   normF,
                                   maxIters)

    # Create the stepper
    stepper = LOCA.Stepper(grp, comboOR, paramList)
    status = stepper.run()
    grp = stepper.getSolutionGroup()
    finalSolution = grp.getX()

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(ChanTestCase))

    # Run the test suite
    if (verbose):
        print "\n***** Running Tests *****\n"
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        result = unittest.TextTestRunner(verbosity=0).run(suite)

    if result.wasSuccessful():
        return 0
    else:
        return 1

if (__name__ == "__main__"):

    # Run the tests
    runTests()




