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
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;

  alpha = alpha / scale;

  try {

    // Create output file to save solutions
    ofstream outFile("ChanContinuation.dat");
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

    // Create parameter list
    NOX::Parameter::List paramList;

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.setParameter("Continuation Method", "Natural");
    stepperList.setParameter("Continuation Method", "Arc Length");
    stepperList.setParameter("Continuation Parameter", "alpha");
    stepperList.setParameter("Initial Value", alpha);
    stepperList.setParameter("Max Value", 5.0/scale);
    stepperList.setParameter("Min Value", 0.0/scale);
    stepperList.setParameter("Max Steps", 50);
    stepperList.setParameter("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.setParameter("Enable Arc Length Scaling", true);
    stepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
    stepperList.setParameter("Initial Scale Factor", 1.0);
    stepperList.setParameter("Min Scale Factor", 1.0e-8);
    stepperList.setParameter("Enable Tangent Factor Step Size Scaling",true);
    stepperList.setParameter("Min Tangent Factor", -1.0);
    stepperList.setParameter("Tangent Factor Exponent",1.0);
    stepperList.setParameter("Compute Eigenvalues",false);

    // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "None");

    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    predictorList.setParameter("Method", "Tangent");
    //predictorList.setParameter("Method", "Secant");

    NOX::Parameter::List& firstStepPredictorList = 
      predictorList.sublist("First Step Predictor");
    firstStepPredictorList.setParameter("Method", "Constant");

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.setParameter("Method", "Constant");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.1/scale);
    stepSizeList.setParameter("Min Step Size", 1.0e-3/scale);
    stepSizeList.setParameter("Max Step Size", 10.0/scale);
    stepSizeList.setParameter("Aggressiveness", 0.5); // for adaptive
    stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
    stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    locaUtilsList.setParameter("Output Information", 
			       LOCA::Utils::Warning +
			       LOCA::Utils::StepperIteration +
   			       LOCA::Utils::StepperDetails +
			       LOCA::Utils::Solver +
			       LOCA::Utils::Parameters +
			       LOCA::Utils::SolverDetails);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			  NOX::Utils::Details +
			  NOX::Utils::OuterIteration + 
			  NOX::Utils::InnerIteration + 
			  NOX::Utils::Warning);

    // Set up the status tests
    NOX::StatusTest::NormF normF(1.0e-8);
    NOX::StatusTest::MaxIters maxIters(maxNewtonIters);
    NOX::StatusTest::Combo comboOR(NOX::StatusTest::Combo::OR, 
				   normF, 
				   maxIters);

    // Create the stepper  
    LOCA::NewStepper stepper(grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished)
      cout << "Stepper failed to converge!" << endl;

    // Get the final solution from the stepper
    const LOCA::LAPACK::Group& finalGroup = 
      dynamic_cast<const LOCA::LAPACK::Group&>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution = 
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup.getX());

    // Output the parameter list
    if (LOCA::Utils::doPrint(LOCA::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      stepper.getParameterList().print(cout);
      cout << endl;
    }

    outFile.close();
  }

  catch (string& s) {
    cout << s << endl;
  }
  catch (char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
