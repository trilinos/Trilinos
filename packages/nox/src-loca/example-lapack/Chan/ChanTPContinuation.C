// $Id$ 
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ChanProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 4.0;
    double beta = 0.0;
    double scale = 1.0;
    int maxNewtonIters = 20;
    //int maxNewtonIters = 6;

    // Create output file to save solutions
    ofstream outFile("ChanTPContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;

    // Set up the problem interface
    ChanProblemInterface chan(n, alpha, beta, scale, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
     LOCA::LAPACK::Group grp(chan);
    grp.setParams(p);

    // Create initial guess for the null vector of jacobian
    NOX::LAPACK::Vector nullVec(n);  // length n
    nullVec.init(1.0);             // initial value 1.0

    // Create parameter list
    NOX::Parameter::List paramList;

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.setParameter("Continuation Method", "Natural");
    stepperList.setParameter("Continuation Method", "Arc Length");
    stepperList.setParameter("Continuation Parameter", "beta");
    stepperList.setParameter("Initial Value", beta);
    stepperList.setParameter("Max Value", 1.0);
    stepperList.setParameter("Min Value", 0.0);
    stepperList.setParameter("Max Steps", 20);
    stepperList.setParameter("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.setParameter("Enable Arc Length Scaling", true);
    stepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
    stepperList.setParameter("Initial Scale Factor", 1.0);
    stepperList.setParameter("Min Scale Factor", 1.0e-8);
    stepperList.setParameter("Enable Tangent Factor Step Size Scaling",false);
    stepperList.setParameter("Min Tangent Factor", -1.0);
    stepperList.setParameter("Tangent Factor Exponent",1.0);

    // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "Turning Point");
    bifurcationList.setParameter("Bifurcation Parameter", "alpha");
    bifurcationList.setParameter("Length Normalization Vector", 
			 dynamic_cast<NOX::Abstract::Vector*>(&nullVec));
    bifurcationList.setParameter("Initial Null Vector",
			 dynamic_cast<NOX::Abstract::Vector*>(&nullVec));

    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    predictorList.setParameter("Method", "Secant");
    //predictorList.setParameter("Method", "Random");
    //predictorList.setParameter("Epsilon", 1.0e-3);

    NOX::Parameter::List& firstStepPredictor 
      = predictorList.sublist("First Step Predictor");
    firstStepPredictor.setParameter("Method", "Random");
    firstStepPredictor.setParameter("Epsilon", 1.0e-3);

    NOX::Parameter::List& lastStepPredictor 
      = predictorList.sublist("Last Step Predictor");
    lastStepPredictor.setParameter("Method", "Random");
    lastStepPredictor.setParameter("Epsilon", 1.0e-3);

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.1);
    stepSizeList.setParameter("Min Step Size", 1.0e-3);
    stepSizeList.setParameter("Max Step Size", 1.0);
    stepSizeList.setParameter("Aggressiveness", 0.5);
    stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
    stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    locaUtilsList.setParameter("Output Information", 
			       LOCA::Utils::Warning +
			       LOCA::Utils::StepperIteration +
			       LOCA::Utils::StepperDetails +
			       LOCA::Utils::Solver +
			       LOCA::Utils::SolverDetails);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			       NOX::Utils::OuterIteration + 
			       //NOX::Utils::OuterIterationStatusTest + 
			       //NOX::Utils::InnerIteration +
			       //NOX::Utils::Parameters +
			       //NOX::Utils::Details + 
			       NOX::Utils::Warning);

    // Create the "Line Search" sublist for the "Line Search Based" solver
    NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
    searchParams.setParameter("Method", "Full Step");

    // Create the newton and  linear solver parameters sublist
  NOX::Parameter::List& directionParameters = nlParams.sublist("Direction");
  NOX::Parameter::List& newtonParameters = directionParameters.sublist("Newton");
  NOX::Parameter::List& linearSolverParameters = newtonParameters.sublist("Linear Solver");

    // Set up the status tests
    NOX::StatusTest::NormF statusTestA(1.0e-5, NOX::StatusTest::NormF::Scaled);
    NOX::StatusTest::MaxIters statusTestB(maxNewtonIters);
    NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

    // Create the stepper  
    LOCA::Stepper stepper(grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished)
      cout << "Stepper failed to converge!" << endl;

    // Output the parameter list
    if (NOX::Utils::doPrint(NOX::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      stepper.getParameterList().print(cout);
      cout << endl;
    }

    outFile.close();
  }

  catch (char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
