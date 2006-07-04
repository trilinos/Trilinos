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
#include "NOX_TestCompare.H"

int main(int argc, char *argv[])
{
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;
  int ierr = 0;

  alpha = alpha / scale;

  try {

    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Method", "Arc Length");
    stepperList.set("Continuation Parameter", "alpha");
    stepperList.set("Initial Value", alpha);
    stepperList.set("Max Value", 5.0/scale);
    stepperList.set("Min Value", 0.0/scale);
    stepperList.set("Max Steps", 50);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.set("Enable Arc Length Scaling", true);
    stepperList.set("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.set("Max Arc Length Parameter Contribution", 0.7);
    stepperList.set("Initial Scale Factor", 1.0);
    stepperList.set("Min Scale Factor", 1.0e-8);
    stepperList.set("Enable Tangent Factor Step Size Scaling",true);
    stepperList.set("Min Tangent Factor", -1.0);
    stepperList.set("Tangent Factor Exponent",1.0);
    stepperList.set("Compute Eigenvalues",false);

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.set("Method", "Constant");
    predictorList.set("Method", "Tangent");
    //predictorList.set("Method", "Secant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.set("Method", "Constant");
    stepSizeList.set("Method", "Adaptive");
    stepSizeList.set("Initial Step Size", 0.1/scale);
    stepSizeList.set("Min Step Size", 1.0e-3/scale);
    stepSizeList.set("Max Step Size", 10.0/scale);
    stepSizeList.set("Aggressiveness", 0.5); // for adaptive
    stepSizeList.set("Failed Step Reduction Factor", 0.5);
    stepSizeList.set("Successful Step Increase Factor", 1.26); // for constant

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    nlParams.set("Nonlinear Solver", "Line Search Based");

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    if (verbose)
       nlPrintParams.set("Output Information", 
				  NOX::Utils::Error +
				  NOX::Utils::Details +
				  NOX::Utils::OuterIteration + 
				  NOX::Utils::InnerIteration + 
				  NOX::Utils::Warning +
				  NOX::Utils::TestDetails + 
				  NOX::Utils::StepperIteration +
				  NOX::Utils::StepperDetails);
     else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp =
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));
    
    grp->setParams(p);

    // Set up the status tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> normF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxIters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Generic> comboOR =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
					      normF,
					      maxIters));

    // Create the stepper  
    LOCA::NewStepper stepper(globalData, grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    // Check for convergence
    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr = 1;
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RefCountPtr<const LOCA::LAPACK::Group> finalGroup =
      Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution = 
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(), 
				 *(globalData->locaUtils));

    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() 
	<< std::endl 
	<< "***** Checking solution statistics *****" 
	<< std::endl;
  
    // Check number of steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 32;
    ierr += testCompare.testValue(numSteps, numSteps_expected, 0.0,
				  "number of continuation steps",
				  NOX::TestCompare::Absolute);

    // Check number of failed steps
    int numFailedSteps = stepper.getNumFailedSteps();
    int numFailedSteps_expected = 0;
    ierr += testCompare.testValue(numFailedSteps, numFailedSteps_expected, 0.0,
				  "number of failed continuation steps",
				  NOX::TestCompare::Absolute);

    // Check final value of continuation parameter
    double alpha_final = finalGroup->getParam("alpha");
    double alpha_expected = 5.0;
    ierr += testCompare.testValue(alpha_final, alpha_expected, 1.0e-14,
				  "final value of continuation parameter", 
				  NOX::TestCompare::Relative);
 
    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected = 203.1991024;
    ierr += testCompare.testValue(norm_x, norm_x_expected, 1.0e-7,
				  "norm of final solution",
				  NOX::TestCompare::Relative);

    destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
    ierr = 1;
  }
  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

  if (ierr == 0)
    cout << "All tests passed!" << endl;
  else
    cout << ierr << " test(s) failed!" << endl;

  return ierr;
}
