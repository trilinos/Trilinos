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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

/*  Notes:
**
**  (.) The TensorBased solver is virtually identical to the
**  LineSearchBased solver.  Thus, maybe at some point I should remove
**  this TensorBased solver and convert LineSearchBased solver and the
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

#include "NOX_Solver_TensorBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member

#include "stdio.h"  // for printf()


NOX::Solver::TensorBased::TensorBased(NOX::Abstract::Group& xgrp,
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
  direction(utils, paramsPtr->sublist("Direction"))    // initialize direction
{
  init();
}

// Protected
void NOX::Solver::TensorBased::init()
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
    cout << "NOX::Solver::TensorBased::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testptr->checkStatus(*this);
  if ((status == NOX::StatusTest::Converged) &&
      (utils.isPrintProcessAndType(NOX::Utils::Warning)))  {
    cout << "Warning: NOX::Solver::TensorBased::init() - The solution passed "
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


bool NOX::Solver::TensorBased::reset(NOX::Abstract::Group& xgrp,
				     NOX::StatusTest::Generic& t,
				     NOX::Parameter::List& p)
{
  solnptr = &xgrp;
  testptr = &t;
  paramsPtr = &p;
  utils.reset(paramsPtr->sublist("Printing"));
  lineSearch.reset(paramsPtr->sublist("Line Search"));
  direction.reset(paramsPtr->sublist("Direction"));
  init();

  return true;
}

bool NOX::Solver::TensorBased::reset(NOX::Abstract::Group& xgrp,
				     NOX::StatusTest::Generic& t)
{
  solnptr = &xgrp;
  testptr = &t;
  init();
  return true;
}

NOX::Solver::TensorBased::~TensorBased() 
{
  delete oldsolnptr;
  delete dirptr;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType  NOX::Solver::TensorBased::iterate()
{
  // First check status
  if (status != NOX::StatusTest::Unconverged) 
    return status;

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnptr;
  NOX::StatusTest::Generic& test = *testptr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction.compute(dir, soln, *this);
  if (!ok) {
    cout << "NOX::Solver::TensorBased::iterate - "
	 << "unable to calculate direction" << endl;
    status = NOX::StatusTest::Failed;
    return status;
  }

  // Update iteration count.
  niter ++;

  // Copy current soln to the old soln.
  oldsoln = soln;

  // Do line search and compute new soln.
  ok = lineSearch.compute2(soln, step, dir, *this, direction);
  if (!ok) {
    if (step == 0) {
      cout << "NOX::Solver::TensorBased::iterate - line search failed" << endl;
      status = NOX::StatusTest::Failed;
      return status;
    }
    else if (utils.isPrintProcessAndType(NOX::Utils::Warning))
       cout << "NOX::Solver::TensorBased::iterate - "
	   << "using recovery step for line search" << endl;
  }
      

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)  {
    cout << "NOX::Solver::LineSearchBased::iterate - "
	 << "unable to compute F" << endl;
    status = NOX::StatusTest::Failed;
    return status;
  }

  
  // Evaluate the current status.
  status = test.checkStatus(*this);
 
  // Return status.
  return status;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBased::solve()
{
  printUpdate();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged) {
    status = iterate();
    printUpdate();
  }

  NOX::Parameter::List& outputParams = paramsPtr->sublist("Output");
  outputParams.setParameter("Nonlinear Iterations", niter);
  outputParams.setParameter("2-Norm of Residual", solnptr->getNormF());

  return status;
}

const NOX::Abstract::Group&
NOX::Solver::TensorBased::getSolutionGroup() const
{
  return *solnptr;
}

const NOX::Abstract::Group&
NOX::Solver::TensorBased::getPreviousSolutionGroup() const
{
  return oldsoln;
}

int NOX::Solver::TensorBased::getNumIterations() const
{
  return niter;
}

const NOX::Parameter::List&
NOX::Solver::TensorBased::getParameterList() const
{
  return *paramsPtr;
}


// protected
void NOX::Solver::TensorBased::printUpdate() 
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
  paramsPtr(NULL),
  print(u)
{
  reset(params);
}

NOX::LineSearch::Tensor::~Tensor()
{

}

bool NOX::LineSearch::Tensor::reset(Parameter::List& params)
{ 
  NOX::Parameter::List& lsparams = params.sublist("Tensor");
  
  // Initialize linesearch parameters for this object
  minStep = lsparams.getParameter("Minimum Step", 1.0e-12);
  defaultStep = lsparams.getParameter("Default Step", 1.0);
  recoveryStep = lsparams.getParameter("Recovery Step", 0.0); // force exit on linesearch failure
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
    if (print.isPrintProcessAndType(NOX::Utils::Warning)) {
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
    cout << "Warning: NOX::Solver::TensorBased::init() - the choice of "
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

  counter.reset();

  return true;
}


bool NOX::LineSearch::Tensor::compute(NOX::Abstract::Group& newGrp,
				      double& step, 
				      const NOX::Abstract::Vector& dir,
				      const NOX::Solver::Generic& s)
{
  return false;
}


bool NOX::LineSearch::Tensor::compute2(NOX::Abstract::Group& newGrp,
				       double& step, 
				       const NOX::Abstract::Vector& dir,
				       const NOX::Solver::Generic& s,
				       NOX::Direction::Tensor& direction) 
{
  bool ok;
  counter.incrementNumLineSearches();
  isNewtonDirection = false;
  
  if (counter.getNumLineSearches() == 1  ||  lsType == Newton)
    isNewtonDirection = true;

  // Do line search and compute new soln.
  if (lsType != Dual || isNewtonDirection)
    ok = performLinesearch(newGrp, step, dir, s, direction);
  else if (lsType == Dual) {
    const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
    double fprime = slopeObj.computeSlope(dir, oldGrp);
    double tensorf;
    double tensorStep;
    bool isTensorDescent = false;
    
    if (fprime < 0) {
      ok = performLinesearch(newGrp, step, dir, s, direction);
      tensorf = 0.5*newGrp.getNormF()*newGrp.getNormF();
      tensorStep = step;
      isTensorDescent = true;
    }      
    const NOX::Abstract::Vector& dir2 = direction.getNewton();

    ok = performLinesearch(newGrp, step, dir2, s, direction);
    double newf = 0.5*newGrp.getNormF()*newGrp.getNormF();

    // If backtracking on the tensor step produced a better step, then use it.
    if (isTensorDescent  &&  tensorf < newf) {
      newGrp.computeX(oldGrp, dir, tensorStep);
      newGrp.computeF();    
    }
  }
#ifdef DEPRECATED_CODE
  else if (lsType == Dual) 
    ok = performDualLinesearch(newGrp, step, dir, s, direction);
#endif
  
  return ok;
}


bool NOX::LineSearch::Tensor::performLinesearch(NOX::Abstract::Group& newsoln,
					double& step,
					const NOX::Abstract::Vector& dir,
					const NOX::Solver::Generic& s,
					NOX::Direction::Tensor& direction)
{
  if (print.isPrintProcessAndType(NOX::Utils::InnerIteration)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Tensor Line Search ("
	 << paramsPtr->sublist("Tensor").getParameter("Submethod","Curvilinear")
	 << ") -- \n";
  }

  // Local variables
  bool isFailed = false;
  bool isAcceptable = false;
  bool isFirstPass = true;
  string message = "(STEP ACCEPTED!)";

  // Set counters
  int lsIterations = 1;

  // Get Old F
  const Abstract::Group& oldsoln = s.getPreviousSolutionGroup();
  double oldf = 0.5*oldsoln.getNormF()*oldsoln.getNormF();  

  // Compute first trial point and its function value
  step = defaultStep;
  newsoln.computeX(oldsoln, dir, step);
  newsoln.computeF();    
  double newf = 0.5*newsoln.getNormF()*newsoln.getNormF();  

  // Compute directional derivative
  double fprime;
  if (lsType == Curvilinear)
    fprime = slopeObj.computeSlope(direction.getNewton(), oldsoln);
  else
    fprime = slopeObj.computeSlope(dir, oldsoln);

  // Compute the convergence criteria for the line search 
  double threshold = oldf + alpha*step*fprime;
  isAcceptable = (newf < threshold);
    
  // Update counter if a linesearch is needed
  if (!isAcceptable) 
    counter.incrementNumNonTrivialLineSearches();

  // Iterate until the trial point is accepted....
  while (!isAcceptable) {  
      
    // Check for linesearch failure
    if (lsIterations > maxIters) {
      isFailed = true;
      message = "(FAILED - Max Iters)";
      break;
    }

    print.printStep(lsIterations, step, oldf, newf);

    // Is the full tensor step a descent direction?  If not, switch to Newton
    if ((!isNewtonDirection) && (isFirstPass && fprime >= 0)) {
      // dir = oldsoln.getNewton();   // bwb - for when Newton put in group
      fprime = slopeObj.computeSlope(direction.getNewton(), oldsoln);
      
      printf("  Switching to Newton.  New fprime = %e\n", fprime);
    }
    else {
      step = selectLambda(newf, oldf, fprime, step);
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
      // bwb - move creation of dir2 outside so only done once.
      NOX::Abstract::Vector* dir2 = dir.clone(ShapeCopy);
      bool ok = direction.computeCurvilinearStep(*dir2, newsoln, s, step);
      newsoln.computeX(oldsoln, *dir2, 1.0);
      delete dir2;
    }
    else {
      newsoln.computeX(oldsoln, dir, step);
    }
    newsoln.computeF();    
    newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

    // Recompute convergence criteria based on new step
    threshold = oldf + alpha*step*fprime;
    isAcceptable = (newf < threshold);
  }


  if (isFailed) {
    counter.incrementNumFailedLineSearches();
    step = recoveryStep;

    if (step != 0.0) {
      // Update the group using Newton direction and recovery step
      newsoln.computeX(oldsoln, direction.getNewton(), step);
      newsoln.computeF();    
      newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

      message = "(USING RECOVERY STEP!)";
    }
  }
  
  print.printStep(lsIterations, step, oldf, newf, message);
  counter.setValues(*paramsPtr);

  return (!isFailed);
}


#ifdef DEPRECATED_CODE
bool NOX::LineSearch::Tensor::performDualLinesearch(NOX::Abstract::Group& newsoln,
						    double& step,
						    const NOX::Abstract::Vector& dir,
						    const NOX::Solver::Generic& s,
						    NOX::Direction::Tensor& direction)
{
  if (print.isPrintProcessAndType(NOX::Utils::InnerIteration)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Dual Direction Linesearch -- \n";
  }
  
  // Set counter to 1
  int lsIterations = 1;

  // local variables
  bool ok;
  bool isFailed = false;
  bool lsRequired = false;
  double tensorf = 1e29;
  double tensors = defaultStep;

  // Get Old F
  const Abstract::Group& oldsoln = s.getPreviousSolutionGroup();
  double oldf = 0.5*oldsoln.getNormF()*oldsoln.getNormF();  

  // Compute new trial point and its function value
  step = defaultStep;
  newsoln.computeX(oldsoln, dir, step);
  newsoln.computeF();    
  double newf = 0.5*newsoln.getNormF()*newsoln.getNormF();  

  // Computation of directional derivative used in curvature condition
  NOX::Abstract::Vector* tmpvecptr = oldsoln.getX().clone(ShapeCopy);
  oldsoln.applyJacobian(dir, *tmpvecptr);
  double fprime = tmpvecptr->dot(oldsoln.getF());

  double slope = slopeObj.computeSlope(dir, oldsoln);
  printf("slope:   %16e   %16e\n", fprime, slope);

  // Compute the convergence criteria for the line search 
  double convergence = oldf + alpha*step*fprime;

  // Increment the counter of nonlinear iterations requiring a linesearch
  if (newf >= convergence) {
    counter.incrementNumNonTrivialLineSearches();
    lsRequired = true;
  }
    
  // If the tensor step is a descent direction (fprime < 0),
  // then backtrack until Armijo-Goldstein condition is satisfied.
  while (newf >= convergence  &&  fprime < 0) {  

    print.printStep(lsIterations, step, oldf, newf);

    counter.incrementNumIterations();
    lsIterations++;
    
    // Backtrack by reducing the step
    step = selectLambda(newf, oldf, fprime, step);

    // Check for linesearch failure: if true, exit the method
    if ((step < minStep) || (lsIterations > maxIters)) {

      counter.incrementNumFailedLineSearches();
      step = recoveryStep;
      
      // Update the group, and exit linesearch
      newsoln.computeX(oldsoln, dir, step);
      newsoln.computeF();    
      newf = 0.5*newsoln.getNormF()*newsoln.getNormF();
      
      // Print out info
      string message = "(USING RECOVERY STEP!)";
      print.printStep(lsIterations, step, oldf, newf, message);
      
      isFailed = true;
      return (!isFailed);
    }

    // Compute new trial point and its function value
    newsoln.computeX(oldsoln, dir, step);
    newsoln.computeF();    
    newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

    // Recompute convergence criteria based on new step
    convergence = oldf + alpha*step*fprime;

  }  // end while loop


  if (lsRequired) {

    string message = "(STEP ACCEPTED!)";
    print.printStep(lsIterations, step, oldf, newf, message);

    // If tensor step is not descent, then use oldf to force acceptance
    // of Newton backtrack...
    if ( !(fprime < 0) )
      newf = oldf;

    // Save info from tensor backtrack and then try backtracking
    // on the Newton step and compare with tensor backtrack.
    tensorf = newf;
    tensors = step;
  
    // Computation of directional derivative used in curvature condition
    //const NOX::Abstract::Vector& dir2 = oldsoln.getNewton();// bwb - dNewton 
    const NOX::Abstract::Vector& dir2 = direction.getNewton();
    oldsoln.applyJacobian(dir2,*tmpvecptr);
    fprime = tmpvecptr->dot(oldsoln.getF());

    slope = slopeObj.computeSlope(dir2, oldsoln);
    printf("slope:   %16e   %16e\n", fprime, slope);

    // Compute new trial point and its function value
    step = defaultStep;
    newsoln.computeX(oldsoln, dir2, step);
    newsoln.computeF();    
    newf = 0.5*newsoln.getNormF()*newsoln.getNormF();  

    counter.incrementNumIterations();
    lsIterations ++;
  
    // Compute the convergence criteria for the line search 
    convergence = oldf + alpha*step*fprime;

    //Iterate until Armijo-Goldstein condition is satisfied
    while (newf >= convergence) {  
    
      print.printStep(lsIterations, step, oldf, newf);

      counter.incrementNumIterations();
      lsIterations ++;
    
      // Backtrack by reducing the step
      step = selectLambda(newf, oldf, fprime, step);
    
      // Check for linesearch failure: if true, exit the method
      if ((step < minStep) || (lsIterations > maxIters)) {
	
	counter.incrementNumFailedLineSearches();
	step = recoveryStep;
	
	// Update the group, and exit linesearch
	newsoln.computeX(oldsoln, dir2, step);
	newsoln.computeF();    
	newf = 0.5*newsoln.getNormF()*newsoln.getNormF();
	
	// Print out info
	string message = "(USING RECOVERY STEP!)";
	print.printStep(lsIterations, step, oldf, newf, message);

	isFailed = true;
	return (!isFailed);
      }

      // Compute new trial point and its function value
      newsoln.computeX(oldsoln, dir2, step);
      newsoln.computeF();    
      newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

      // Recompute convergence criteria based on new step
      convergence = oldf + alpha*step*fprime;

    }  // end while loop
  } // end if
  
  string message = "(STEP ACCEPTED!)";
  print.printStep(lsIterations, step, oldf, newf, message);


  // If backtracking on the tensor step produced a better step, then use it.
  if (lsRequired  &&  tensorf < newf) {
    newsoln.computeX(oldsoln, dir, tensors);
    newsoln.computeF();    
    newf = 0.5*newsoln.getNormF()*newsoln.getNormF();
  }

  
#if DEBUG_LEVEL > 1
  cout << "Current point: ";
  newsoln.getX().print();
#endif

  ok = !isFailed;

  counter.setValues(*paramsPtr);
  //counter.setValues(params);

  delete tmpvecptr;

  return ok;
}
#endif // DEPRECATED_CODE

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


#endif
