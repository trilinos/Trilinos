// $Id$ 
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
#include "PitchforkProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 1.0;
    double beta = 0.0;
    double lambda = -2.25;
    double pi = 4.0*atan(1.0);
    double h = 2.0/(n-1);
    int maxNewtonIters = 20;

    // Set up the problem interface
    PitchforkProblemInterface pf(n, alpha, beta, lambda);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("lambda",lambda);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    LOCA::LAPACK::Group grp(pf);
    grp.setParams(p);

    // Create asymmetry vector
    NOX::LAPACK::Vector asymVec(n);  // length n
    for (int i=0; i<n; i++)
      asymVec(i) = sin( pi/2.0 * (-1.0 + h*i) );

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
    stepperList.setParameter("Min Value", -0.01);
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
    stepperList.setParameter("Compute Eigenvalues",false);

    // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "Pitchfork");
    bifurcationList.setParameter("Bifurcation Parameter", "lambda");
    bifurcationList.setParameter("Length Normalization Vector", 
			 dynamic_cast<NOX::Abstract::Vector*>(&asymVec));
    bifurcationList.setParameter("Initial Null Vector",
			 dynamic_cast<NOX::Abstract::Vector*>(&asymVec));
    bifurcationList.setParameter("Asymmetric Vector",
			 dynamic_cast<NOX::Abstract::Vector*>(&asymVec));

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
			       //NOX::Utils::OuterIteration + 
			       //NOX::Utils::OuterIterationStatusTest + 
			       //NOX::Utils::InnerIteration +
			       //NOX::Utils::Parameters +
			       //NOX::Utils::Details + 
			       NOX::Utils::Warning);

    // Create the "Line Search" sublist for the "Line Search Based" solver
    NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
    searchParams.setParameter("Method", "Full Step");

    // Direction sublist
    NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
    dirParams.setParameter("Method", "Newton");

    // Newton sublist
    NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");

    // Linear solver sublist
    NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");

    // Set up the status tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> statusTestA = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, 
					      NOX::StatusTest::NormF::Scaled));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> statusTestB = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, statusTestA, 
				 statusTestB);

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
  }

  catch (char* s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
