#! /usr/bin/env python

import setpath

import NOX
import LOCA

def main():

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
    locaUtilsList.setParameter("Output Information",0)

    # Create the "Solver" parameters sublist to be used with NOX Solvers
    nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based")

    nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information",0)

    # Set up the status tests
    normF = NOX.StatusTest.NormF(1.0e-8)
    maxIters = NOX.StatusTest.MaxIters(maxNewtonIters)
    comboOR = NOX.StatusTest.Combo(NOX.StatusTest.Combo.OR, 
                                   normF, 
                                   maxIters)

    # Create the stepper
    stepper = LOCA.Stepper(grp, comboOR, paramList)
    
    # Perform continuation run
    status = stepper.run()
    
    if (status != LOCA.Iterator.Finished):
        return 1
    else:
        return 0

if (__name__ == "__main__"):
    main()
    

