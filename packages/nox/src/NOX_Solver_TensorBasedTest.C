// $Id$ 
// $Source$ 

#ifdef WITH_PRERELEASE
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
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

/*  Notes:
**
**  (.) The TensorBasedTest solver is virtually identical to the
**  LineSearchBased solver.  Thus, maybe at some point I should remove
**  this TensorBasedTest solver and convert LineSearchBased solver and the
**  LineSearch classes to work with the tensor linesearches.  Right
**  now I see 2 options for this conversion: Make optional argument in
**  LineSearch::compute to allow for a direction argument or add
**  getDirection method to solver so that linesearch object can
**  compute the curvilinear linesearch.  The latter option might have
**  trouble using the const direction.  Need to investigate...
**     //NOX::Abstract::Vector dir2 = dir.clone(ShapeCopy);
**     //const NOX::Direction::Tensor& direction = s.getDirection();
**  
**  (.)  Should change to *sufficient* decrease condition instead of
**  just "fprime<0"
**
**  (.)  Maybe move the test of full step into compute instead of in
**  performLinesearch.  However, this might cause trouble with
**  counters and other things.
**
**  (.)  In the dual linesearch, it is checking both full steps and
**  taking the best of either.  This is different from TENSOLVE.
**
**  (.)  Old comment says:
**  "// Note that for Newton direction, fprime = -2.0*oldf"
**  Is this really true?
*/

#include "NOX_Solver_TensorBasedTest.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member

#include "stdio.h"  // for printf()


NOX::Solver::TensorBasedTest::TensorBasedTest(NOX::Abstract::Group& xgrp,
				      NOX::StatusTest::Generic& t,
				      NOX::Parameter::List& p) :
  solnptr(&xgrp),		// pointer to xgrp
  oldsolnptr(xgrp.clone(DeepCopy)), // create via clone
  oldsoln(*oldsolnptr),		// reference to just-created pointer
  dirptr(xgrp.getX().clone(ShapeCopy)), // create via clone 
  dir(*dirptr),			// reference to just-created pointer
  testptr(&t),			// pointer to t
  paramsPtr(&p),			// copy p
  utils(paramsPtr->sublist("Printing")),               // initialize utils
  lineSearch(utils, paramsPtr->sublist("Line Search")),// initialize linesearch
  direction(utils, paramsPtr->sublist("Direction")),   // initialize direction
  prePostOperator(utils, paramsPtr->sublist("Solver Options"))
{
  init();
}

// Protected
void NOX::Solver::TensorBasedTest::init()
{
  // Initialize 
  step = 0;
  niter = 0;
  status = NOX::StatusTest::Unconverged;

  // Print out initialization information
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(cout,5);
    cout << "\n" << NOX::Utils::fill(72) << "\n";
  }

  // Compute F of initial guess
  NOX::Abstract::Group::ReturnType rtype = solnptr->computeF();
  if (rtype != NOX::Abstract::Group::Ok)    {
    cout << "NOX::Solver::TensorBasedTest::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testptr->checkStatus(*this);
  if ((status == NOX::StatusTest::Converged) &&
      (utils.isPrintProcessAndType(NOX::Utils::Warning)))  {
    cout << "Warning: NOX::Solver::TensorBasedTest::init() - The solution passed "
	 << "into the solver (either through constructor or reset method) "
	 << "is already converged!  The solver will not "
	 << "attempt to solve this system since status is flagged as "
	 << "converged." << endl;
  }

  // Print out status tests
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters))  {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testptr->print(cout, 5);
    cout <<"\n" << NOX::Utils::fill(72) << "\n";
  }
}


bool NOX::Solver::TensorBasedTest::reset(NOX::Abstract::Group& xgrp,
				     NOX::StatusTest::Generic& t,
				     NOX::Parameter::List& p)
{
  solnptr = &xgrp;
  testptr = &t;
  paramsPtr = &p;
  utils.reset(paramsPtr->sublist("Printing"));
  lineSearch.reset(paramsPtr->sublist("Line Search"));
  direction.reset(paramsPtr->sublist("Direction"));
  prePostOperator.reset(utils, paramsPtr->sublist("Solver Options"));
  init();

  return true;
}

bool NOX::Solver::TensorBasedTest::reset(NOX::Abstract::Group& xgrp,
				     NOX::StatusTest::Generic& t)
{
  solnptr = &xgrp;
  testptr = &t;
  init();
  return true;
}

NOX::Solver::TensorBasedTest::~TensorBasedTest() 
{
  delete oldsolnptr;
  delete dirptr;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBasedTest::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType  NOX::Solver::TensorBasedTest::iterate()
{
  prePostOperator.runPreIterate(*this);

  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    prePostOperator.runPostIterate(*this);
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnptr;
  NOX::StatusTest::Generic& test = *testptr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction.compute(dir, soln, *this);
  if (!ok) {
    cout << "NOX::Solver::TensorBasedTest::iterate - "
	 << "unable to calculate direction" << endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }

  // Update iteration count.
  niter ++;

  // Copy current soln to the old soln.
  oldsoln = soln;

  // Do line search and compute new soln.
  //ok = lineSearch.compute2(soln, step, dir, *this, direction);
  ok = lineSearch.compute(soln, step, dir, *this);
  if (!ok) {
    if (step == 0) {
      cout << "NOX::Solver::TensorBasedTest::iterate - line search failed" << endl;
      status = NOX::StatusTest::Failed;
      prePostOperator.runPostIterate(*this);
      return status;
    }
    else if (utils.isPrintProcessAndType(NOX::Utils::Warning))
       cout << "NOX::Solver::TensorBasedTest::iterate - "
	    << "using recovery step for line search" << endl;
  }
      

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)  {
    cout << "NOX::Solver::LineSearchBased::iterate - "
	 << "unable to compute F" << endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }

  
  // Evaluate the current status.
  status = test.checkStatus(*this);
 
  prePostOperator.runPostIterate(*this);

  // Return status.
  return status;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBasedTest::solve()
{
  prePostOperator.runPreSolve(*this);

  printUpdate();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged) {
    status = iterate();
    printUpdate();
  }

  NOX::Parameter::List& outputParams = paramsPtr->sublist("Output");
  outputParams.setParameter("Nonlinear Iterations", niter);
  outputParams.setParameter("2-Norm of Residual", solnptr->getNormF());

  prePostOperator.runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group&
NOX::Solver::TensorBasedTest::getSolutionGroup() const
{
  return *solnptr;
}

const NOX::Abstract::Group&
NOX::Solver::TensorBasedTest::getPreviousSolutionGroup() const
{
  return oldsoln;
}

int NOX::Solver::TensorBasedTest::getNumIterations() const
{
  return niter;
}

const NOX::Parameter::List&
NOX::Solver::TensorBasedTest::getParameterList() const
{
  return *paramsPtr;
}

const NOX::Direction::Tensor&
NOX::Solver::TensorBasedTest::getDirection() const
{
  return direction;
}

// protected
void NOX::Solver::TensorBasedTest::printUpdate() 
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested  
  if ((status == NOX::StatusTest::Unconverged) &&
      (utils.isPrintProcessAndType(NOX::Utils::OuterIterationStatusTest))) {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Status Test Results --\n";    
    testptr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utils.isPrintType(NOX::Utils::InnerIteration)) {
    normSoln = solnptr->getNormF();
    normStep = (niter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utils.isPrintProcessAndType(NOX::Utils::OuterIteration)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Nonlinear Solver Step " << niter << " -- \n";
    cout << "f = " << utils.sciformat(normSoln);
    cout << "  step = " << utils.sciformat(step);
    cout << "  dx = " << utils.sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      cout << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) && 
      (utils.isPrintProcessAndType(NOX::Utils::OuterIteration))) {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testptr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }
}


// ================================================================
//                     NOX::LineSearch::Tensor
// ================================================================


NOX::LineSearch::Tensor::Tensor(const NOX::Utils& u, Parameter::List& params) :
  paramsPtr(&params),
  utils(u),   //paramsPtr->sublist("Printing")),
  print(utils)
{
  //  reset(paramsPtr->sublist("Line Search"));
  reset(*paramsPtr);
}

NOX::LineSearch::Tensor::~Tensor()
{
  printf("multsJv = %d   (linesearch)\n", multsJv);
}

bool NOX::LineSearch::Tensor::reset(NOX::Parameter::List& lsParams)
{
  multsJv = 0;
  
  //NOX::Parameter::List& lsParams = paramsPtr->sublist("Line Search");

  // Determine the specific type of tensor linesearch to perform
  string choice = lsParams.getParameter("Method", "Curvilinear");

  cout << choice << endl;
  
  if (choice == "Curvilinear")
    lsType = Curvilinear;
  else if (choice == "Dual")
    lsType = Dual;
  else if (choice == "Standard")
    lsType = Standard;
  //  else if (choice == "Full Step")
  //    lsType = FullStep;
  else if (choice == "Newton")
    lsType = Newton;
  else
  {
    if (utils.isPrintProcessAndType(NOX::Utils::Error))
      cerr << "NOX::Direction::Tensor::reset() - The choice of "
	   << "\"Line Search\" parameter " << choice
	   << " is invalid." << endl;
    throw "NOX Error";
  }
  //  Copy Method into "Submethod" (temporary hack for data scripts)
  lsParams.setParameter("Submethod", choice);

  // Make a reference to the sublist holding the global strategy parameters
  NOX::Parameter::List& gsParams = lsParams.sublist(choice);

#ifdef CODE_FROM_TENSORBASED  
  // Decide what step to use in case of linesearch failure
  choice = gsParams.getParameter("Recovery Step Type", "Constant");
  if (choice == "Constant")
    recoveryStepType = Constant;          // Use value in "Recovery Step"
  else if (choice == "Last Computed Step") 
    recoveryStepType = LastComputedStep;  // Use last step from linesearch
  else
  {
    cerr << "NOX::Solver::TensorBased::reset() - "
	 << "Invalid \"Recovery Step Type\"" << endl;
    throw "NOX Error";
  }
#endif

  // Initialize linesearch parameters for this object
  minStep = gsParams.getParameter("Minimum Step", 1.0e-12);
  defaultStep = gsParams.getParameter("Default Step", 1.0);
  recoveryStep = gsParams.getParameter("Recovery Step", 0.0); // exit on fail
  maxIters = gsParams.getParameter("Max Iters", 40);
  alpha = gsParams.getParameter("Alpha Factor", 1.0e-4);

  choice = gsParams.getParameter("Lambda Selection", "Halving");
  if (choice == "Halving")
    lambdaSelection = Halving;
  else if (choice == "Quadratic") 
    lambdaSelection = Quadratic;
  else
  {
    if (utils.isPrintProcessAndType(NOX::Utils::Error))
      cerr << "NOX::Solver::TensorBased::reset() - The choice of "
	   << "\"Lambda Selection\" parameter " << choice
	   << " is invalid." << endl;
    throw "NOX Error";
  }

  choice = gsParams.getParameter("Sufficient Decrease Condition",
				 "Armijo-Goldstein");
  if (choice == "Armijo-Goldstein") 
    suffDecrCond = ArmijoGoldstein;     // This is the only one implemented
  else if (choice == "Ared/Pred") 
    suffDecrCond = AredPred;
  else if (choice == "None")
    suffDecrCond = None;
  else
  {
    if (utils.isPrintProcessAndType(NOX::Utils::Error))
      cerr << "NOX::Solver::TensorBased::reset() - The choice of "
	   << "\"Sufficient Decrease Condition\" parameter " << choice
	   << " is invalid." << endl;
    throw "NOX Error";
  }


#ifdef OLD_CODE
  // Initialize linesearch parameters for this object
  minStep = lsparams.getParameter("Minimum Step", 1.0e-12);
  defaultStep = lsparams.getParameter("Default Step", 1.0);
  recoveryStep = lsparams.getParameter("Recovery Step", 0.0); // exit on fail
  maxIters = lsparams.getParameter("Max Iters", 40);
  alpha = lsparams.getParameter("Alpha Factor", 1.0e-4);
  paramsPtr = &params;

  // Do line search and compute new soln.
  string choice = lsparams.getParameter("Submethod", "Curvilinear");

  if (choice == "Curvilinear")
    lsType = Curvilinear;
  else if (choice == "Dual")
    lsType = Dual;
  else if (choice == "Standard")
    lsType = Standard;
  else if (choice == "Newton")
    lsType = Newton;
  else {
    if (utils.isPrintProcessAndType(NOX::Utils::Warning)) {
      cout << "Warning: NOX::Direction::Tensor::reset() - the choice of "
	   << "\"Line Search\" \nparameter is invalid.  Using curvilinear "
	   << "line search." << endl;
    }
    lsparams.setParameter("Submethod", "Curvilinear");
    lsType = Curvilinear;
  }


  choice = lsparams.getParameter("Lambda Selection", "Halving");
  if (choice == "Halving") {
    lambdaSelection = Halving;
  }
  else if (choice == "Quadratic") {
    lambdaSelection = Quadratic;
  }
  else {
    cout << "Warning: NOX::Solver::TensorBasedTest::init() - the choice of "
	 << "\"Lambda Selection\" parameter is invalid." << endl;
    lambdaSelection = Halving;
  }


  choice = lsparams.getParameter("Sufficient Decrease Condition",
				 "Armijo-Goldstein");
  if (choice == "Ared/Pred") 
    convCriteria = AredPred;
  else if (choice == "None")
    convCriteria = None;
  else 
    convCriteria = ArmijoGoldstein;     // bwb - the others aren't implemented
#endif 
  
  counter.reset();

  return true;
}


bool NOX::LineSearch::Tensor::compute(NOX::Abstract::Group& newGrp,
				      double& step, 
				      const NOX::Abstract::Vector& dir,
				      const NOX::Solver::Generic& s)
{
  bool ok;
  counter.incrementNumLineSearches();
  isNewtonDirection = false;

  const NOX::Direction::Tensor& direction =
    (dynamic_cast<const Solver::TensorBasedTest*>(&s))->getDirection();

#ifdef TRIAL_CODE
  // New code added
  const Solver::Generic* test = 0;
  test = dynamic_cast<const Solver::TensorBasedTest*>(&s);
  if (test == 0)
    {
      //printf("Not a TensorBasedTest solver...\n");
    }
  else
    {
      //printf("IS a TensorBasedTest solver...\n");
      direction = (dynamic_cast<const Solver::TensorBasedTest*>(&s))->getDirection();
    }
#endif
  
  
  if (counter.getNumLineSearches() == 1  ||  lsType == Newton)
    isNewtonDirection = true;

  // Do line search and compute new soln.
  if (lsType != Dual || isNewtonDirection)
    ok = performLinesearch(newGrp, step, dir, s, direction);
  else if (lsType == Dual) {
    const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
    double fprime = slopeObj.computeSlope(dir, oldGrp);
    double tensorf = 0.0;
    double tensorStep = 1.0;
    bool isTensorDescent = false;
    
    if (fprime < 0) {
      ok = performLinesearch(newGrp, step, dir, s, direction);
      tensorf = 0.5*newGrp.getNormF()*newGrp.getNormF();
      tensorStep = step;
      isTensorDescent = true;
    }

    //const NOX::Abstract::Vector& dir2 = direction.getNewton();
    //ok = performLinesearch(newGrp, step, dir2, s, direction);
    ok = performLinesearch(newGrp, step, direction.getNewton(), s, direction);
    double newValue = 0.5*newGrp.getNormF()*newGrp.getNormF();

    // If backtracking on the tensor step produced a better step, then use it.
    if (isTensorDescent  &&  tensorf < newValue) {
      newGrp.computeX(oldGrp, dir, tensorStep);
      newGrp.computeF();    
    }
  }
  
  return ok;
}


bool NOX::LineSearch::Tensor::performLinesearch(NOX::Abstract::Group& newsoln,
				double& step,
				const NOX::Abstract::Vector& dir,
				const NOX::Solver::Generic& s,
				const NOX::Direction::Tensor& direction)
{
  if (utils.isPrintProcessAndType(NOX::Utils::InnerIteration)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Tensor Line Search ("
	 << paramsPtr->getParameter("Submethod","Curvilinear")
	 << ") -- \n";
  }

  // Local variables
  NOX::Abstract::Vector* dir2 = NULL;
  bool isFailed = false;
  bool isAccepted = false;
  bool isFirstPass = true;
  string message = "(STEP ACCEPTED!)";

  // Set counters
  int lsIterations = 1;

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  string dirString = s.getParameterList().sublist("Direction").
    getParameter("Method", "Tensor");
  double eta = (suffDecrCond == AredPred) ? 
    s.getParameterList().sublist("Direction").sublist(dirString).
    sublist("Linear Solver").getParameter("Tolerance", -1.0) : 0.0;

  // Get Old function value
  const Abstract::Group& oldsoln = s.getPreviousSolutionGroup();
  double oldValue = 0.5*oldsoln.getNormF()*oldsoln.getNormF();  

  // Compute directional derivative at old solution
  double fprime = (lsType == Curvilinear) ? 
    slopeObj.computeSlope(direction.getNewton(), oldsoln) :
    slopeObj.computeSlope(dir, oldsoln);
  multsJv++;

  // Compute first trial point and its function value
  step = defaultStep;
  newsoln.computeX(oldsoln, dir, step);
  newsoln.computeF();    
  double newValue = 0.5*newsoln.getNormF()*newsoln.getNormF();  

  // Compute the convergence criteria for the line search 
  //  double threshold = oldValue + alpha*step*fprime;
  //  isAccepted = (newValue < threshold);
//   if (fprime >= 0.0) 
//   {
//     printBadSlopeWarning(fprime);
//     isFailed = true;
//   }
//   else 
    isAccepted = checkConvergence(newValue, oldValue, fprime, step, eta,
				  lsIterations, 0);

  // Update counter and allocate memory for dir2 if a linesearch is needed
  if (!isAccepted) { 
    counter.incrementNumNonTrivialLineSearches();
    dir2 = dir.clone(ShapeCopy);
    *dir2 = dir;
  }

  // Iterate until the trial point is accepted....
  while ((!isAccepted) && (!isFailed)) {  
      
    // Check for linesearch failure
    if (lsIterations > maxIters) {
      isFailed = true;
      message = "(FAILED - Max Iters)";
      break;
    }

    print.printStep(lsIterations, step, oldValue, newValue);

    // Is the full tensor step a descent direction?  If not, switch to Newton
    if ((!isNewtonDirection) && (isFirstPass && fprime >= 0)) {
      // dir = oldsoln.getNewton();   // bwb - for when Newton put in group
      *dir2 = direction.getNewton();
      fprime = slopeObj.computeSlope(*dir2, oldsoln);
      multsJv++;
      
      printf("  Switching to Newton.  New fprime = %e\n", fprime);
    }
    else {
      step = selectLambda(newValue, oldValue, fprime, step);
    }
    
    isFirstPass = false;

    // Check for linesearch failure
    if (step < minStep) {
      isFailed = true;
      message = "(FAILED - Min Step)";
      break;
    }

    // Update the number of linesearch iterations
    counter.incrementNumIterations();
    lsIterations ++;

    // Compute new trial point and its function value
    if (lsType == Curvilinear) {

      // bwb - oldsoln needed for preconditioner, test when right pre available
      //bool ok = direction.computeCurvilinearStep2(*dir2, newsoln, s, step);
      //bool ok = direction.computeCurvilinearStep(*dir2, newsoln, s, step);
      direction.computeCurvilinearStep(*dir2, oldsoln, s, step);
      newsoln.computeX(oldsoln, *dir2, 1.0);
    }
    else {
      newsoln.computeX(oldsoln, *dir2, step);
    }
    newsoln.computeF();    
    newValue = 0.5*newsoln.getNormF()*newsoln.getNormF();

    // Recompute convergence criteria based on new step
    // threshold = oldValue + alpha*step*fprime;
    // isAccepted = (newValue < threshold);
    isAccepted = checkConvergence(newValue, oldValue, fprime, step, eta,
				  lsIterations, 0);
  }

  if (isFailed) {
    counter.incrementNumFailedLineSearches();
    step = recoveryStep;

    if (step != 0.0) {
      // Update the group using Newton direction and recovery step
      newsoln.computeX(oldsoln, direction.getNewton(), step);
      newsoln.computeF();    
      newValue = 0.5*newsoln.getNormF()*newsoln.getNormF();

      message = "(USING RECOVERY STEP!)";
    }
  }
  
  print.printStep(lsIterations, step, oldValue, newValue, message);
  counter.setValues(*paramsPtr);

  if (dir2 != NULL) {
    delete dir2;
    dir2 = NULL;
  }

  if (suffDecrCond == AredPred)
    paramsPtr->setParameter("Adjusted Tolerance", 1.0 - step * (1.0 - eta));

  return (!isFailed);
}


bool NOX::LineSearch::Tensor::checkConvergence(double newValue, double oldValue, 
					       double oldSlope,
					       double step, double eta, 
					       int nIters,
					       int nNonlinearIters) const
{
  /*
  if ((nIters == 1) && (doForceInterpolation))
    return false;

  if ((doAllowIncrease) && (nNonlinearIters <= maxIncreaseIter))
  {
    double relativeIncrease = newValue / oldValue;
    if (relativeIncrease < maxRelativeIncrease)
      return true;
  }
  */
  
  switch (suffDecrCond)
  {

  case ArmijoGoldstein:
    
    return (newValue <= oldValue + alpha * step * oldSlope);
    break;

  case AredPred:
    {
      double newEta = 1.0 - step * (1.0 - eta);
      return (newValue <= oldValue * (1.0 - alpha * (1.0 - newEta)));
      break;
    }
  case None:

    return true;
    break;

  default:

    cerr << "NOX::LineSearch::Tensor::checkConvergence - Unknown convergence criteria" << endl;
    throw "NOX Error";

  }
}

double NOX::LineSearch::Tensor::selectLambda(double newf, double oldf,
					     double oldfprime, double lambda)
{
  double lambdaRet;
  double temp;
  
  if (lambdaSelection == Quadratic) {
    temp = -oldfprime / (2.0*(newf - oldf - oldfprime));
    if (temp < 0.1)
      temp = 0.1;
    lambdaRet = temp *lambda;
  }
  else {
    lambdaRet = 0.5 * lambda;
  }
  return lambdaRet;
}


void NOX::LineSearch::Tensor::printBadSlopeWarning(double slope) const
{
  if (print.isPrintProcessAndType(NOX::Utils::Warning))
    cout << "WARNING: Computed slope is positive (slope = " 
	 << slope
	 << ").\n" << "Using recovery step!" 
	 << endl;
}

#endif
