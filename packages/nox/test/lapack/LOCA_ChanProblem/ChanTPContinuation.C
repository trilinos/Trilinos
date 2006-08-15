// $Id$ 
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ChanProblemInterface.H"
#include "NOX_TestCompare.H"
#include "NOX_Random.H"

int main(int argc, char *argv[])
{
  int ierr = 0;

  try {
    int n = 100;
    double alpha = 4.0;
    double beta = 0.0;
    double scale = 1.0;
    int maxNewtonIters = 10;

    NOX::Random::setSeed(1);

    bool verbose = false;
    // Check for verbose output
    if (argc>1) 
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Create initial guess for the null vector of jacobian
    Teuchos::RefCountPtr<NOX::Abstract::Vector> nullVec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    nullVec->init(1.0);               // initial value 1.0

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 20);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Turning Point");
    bifurcationList.set("Bifurcation Parameter", "alpha");
    bifurcationList.set("Length Normalization Vector", nullVec);
    bifurcationList.set("Initial Null Vector", nullVec);

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    Teuchos::ParameterList& firstStepPredictor 
      = predictorList.sublist("First Step Predictor");
    firstStepPredictor.set("Method", "Random");
    firstStepPredictor.set("Epsilon", 1.0e-3);
    Teuchos::ParameterList& lastStepPredictor 
      = predictorList.sublist("Last Step Predictor");
    lastStepPredictor.set("Method", "Random");
    lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 1.0);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    if (verbose)
      nlPrintParams.set("Output Information", 
			NOX::Utils::Error +
			NOX::Utils::OuterIteration + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Details + 
			NOX::Utils::Warning +
			NOX::Utils::TestDetails + 
			NOX::Utils::StepperIteration +
			NOX::Utils::StepperDetails +
			NOX::Utils::StepperParameters);
    else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create LAPACK Factory
    Teuchos::RefCountPtr<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

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
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> statusTestA =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-5, 
					      NOX::StatusTest::NormF::Scaled));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> statusTestB =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      statusTestA, statusTestB));

    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

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
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
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
    int numSteps_expected = 10;
    ierr += testCompare.testValue(numSteps, numSteps_expected, 0.0,
				  "number of continuation steps", 
				  NOX::TestCompare::Absolute);

    // Check number of failed steps
    int numFailedSteps = stepper.getNumFailedSteps();
    int numFailedSteps_expected = 0;
    ierr += testCompare.testValue(numFailedSteps, numFailedSteps_expected, 
				  0.0, "number of failed continuation steps", 
				  NOX::TestCompare::Absolute);

    // Check final value of continuation parameter
    double beta_final = finalGroup->getParam("beta");
    double beta_expected = 0.0;
    ierr += testCompare.testValue(beta_final, beta_expected, 1.0e-14,
				  "final value of continuation parameter", 
				  NOX::TestCompare::Relative);

    // Check final value of turning point parameter
    double alpha_final = finalGroup->getParam("alpha");
    double alpha_expected = 3.1601952;
    ierr += testCompare.testValue(alpha_final, alpha_expected, 1.0e-5,
				  "final value of turning point parameter", 
				  NOX::TestCompare::Relative);

    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected = 63.1872045;
    ierr += testCompare.testValue(norm_x, norm_x_expected, 1.0e-4,
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
