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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Solver_TensorBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Solver;

/* Some compilers (in particular the SGI and ASCI Red - TFLOP) 
 * fail to find the max and min function.  Therfore we redefine them 
 * here. 
 */ 
#ifdef max
#undef max
#endif

#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif

#define min(a,b) ((a)<(b)) ? (a) : (b);

TensorBased::TensorBased(Abstract::Group& xgrp, StatusTest::Generic& t, const Parameter::List& p) :
  solnptr(&xgrp),		// pointer to xgrp
  oldsolnptr(xgrp.clone(DeepCopy)), // create via clone
  oldsoln(*oldsolnptr),		// reference to just-created pointer
  dirptr(xgrp.getX().clone(ShapeCopy)), // create via clone 
  dir(*dirptr),			// reference to just-created pointer
  testptr(&t),			// pointer to t
  params(p),			// copy p
  direction(params.sublist("Direction")) // initialize direction
{
  init();
}

// Protected
void TensorBased::init()
{
  // Initialize 
  step = 0;
  niter = 0;
  totalNumLSSteps = 0;
  totalNumLSIterations = 0;
  totalNumFailedLineSearches = 0;
  status = StatusTest::Unconverged;

  // Initialize linesearch parameters for this solver object
  Parameter::List& lsparams = params.sublist("Line Search");
  minStep = lsparams.getParameter("Minimum Step", 1.0e-12);
  defaultStep = lsparams.getParameter("Default Step", 1.0);
  recoveryStep = lsparams.getParameter("Recovery Step", defaultStep);
  maxIters = lsparams.getParameter("Max Iters", 40);
  alpha = lsparams.getParameter("Alpha Factor", 1.0e-4);

  // Set up utilities (i.e., set print processor, etc)
  Utils::setUtils(params);
  
  // Print out initialization information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Used in Nonlinear Solver --\n\n";
    params.print(cout,5);
    cout << "\n" << Utils::fill(72) << "\n";
  }

  // Compute F of initial guess
  bool ok = solnptr->computeF();
  if (!ok) {
    cout << "NOX::Solver::TensorBased::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testptr->checkStatus(*this);
  if (status == StatusTest::Converged) {
    if (Utils::doPrint(Utils::Warning)) {
      cout << "Warning: NOX::Solver::TensorBased::init() - The solution passed "
	   << "into the solver (either through constructor or reset method) "
	   << "is already converged!  The solver will not "
	   << "attempt to solve this system since status is flagged as "
	   << "converged." << endl;
    }
  }

  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testptr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";
  }

}

bool TensorBased::reset(Abstract::Group& xgrp, StatusTest::Generic& t, const Parameter::List& p) 
{
  solnptr = &xgrp;
  testptr = &t;
  params = p;
  direction.reset(params.sublist("Direction"));
  init();

  return true;
}

TensorBased::~TensorBased() 
{
  delete oldsolnptr;
  delete dirptr;
}


NOX::StatusTest::StatusType TensorBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType TensorBased::iterate()
{
  // First check status
  if (status != StatusTest::Unconverged) 
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnptr;
  StatusTest::Generic& test = *testptr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction.compute(dir, soln, *this);
  if (!ok) {
    cout << "NOX::Solver::TensorBased::iterate - "
	 << "unable to calculate direction" << endl;
    status = StatusTest::Failed;
    return status;
  }

  // Update iteration count.
  niter ++;

  // Copy current soln to the old soln.
  oldsoln = soln;

  // Do line search and compute new soln.
  string choice = params.sublist("Line Search").
    getParameter("Method", "Curvilinear");
  if (choice == "Curvilinear") 
    ok = performCurvilinearLinesearch(); 
  else if (choice == "Standard")
    ok = performStandardLinesearch();
  else {
    cout << "ERROR: NOX::Direction::Tensor::iterate() - the choice of "
	 << "\"Line Search\" parameter is invalid." << endl;
    //throw "NOX Error";
    ok = performCurvilinearLinesearch(); 
  }


  if (!ok) {
    if (step == 0) {
      cout << "NOX::Solver::TensorBased::iterate - linesearch failed" << endl;
      status = StatusTest::Failed;
      return status;
    }
    else if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Solver::TensorBased::iterate - "
	   << "using recovery step for linesearch" << endl;
  }
      

  // Compute F for new current solution.
  ok = soln.computeF();
  if (!ok) {
    cout << "NOX::Solver::TensorBased::iterate - unable to compute F" << endl;
    status = StatusTest::Failed;
    return status;
  }

  // Evaluate the current status.
  status = test.checkStatus(*this);
 
  // Return status.
  return status;
}

NOX::StatusTest::StatusType TensorBased::solve()
{
  printUpdate();

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = iterate();
    printUpdate();
  }

  setOutputParameters();

  return status;
}

const Abstract::Group& TensorBased::getSolutionGroup() const
{
  return *solnptr;
}

const Abstract::Group& TensorBased::getPreviousSolutionGroup() const
{
  return oldsoln;
}

int TensorBased::getNumIterations() const
{
  return niter;
}

const Parameter::List& TensorBased::getParameterList() const
{
  return params;
}


// protected
void TensorBased::printUpdate() 
{
  double norm_soln;
  double norm_step;

  // Print the status test parameters at each iteration if requested  
  if ((status == StatusTest::Unconverged) && 
      (Utils::doPrint(Utils::OuterIterationStatusTest))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Status Test Results --\n";    
    testptr->print(cout);
    cout << Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (Utils::doAllPrint(Utils::OuterIteration)) {
    norm_soln = solnptr->getNormF();
    norm_step = (niter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (Utils::doPrint(Utils::OuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Nonlinear Solver Step " << niter << " -- \n";
    cout << "f = " << Utils::sci(norm_soln);
    cout << "  step = " << Utils::sci(step);
    cout << "  dx = " << Utils::sci(norm_step);
    if (status == StatusTest::Converged)
      cout << " (Converged!)";
    if (status == StatusTest::Failed)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }

  // Print the final parameter values of the status test
  if ((status != StatusTest::Unconverged) && 
      (Utils::doPrint(Utils::OuterIteration))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testptr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}


bool TensorBased::setOutputParameters() 
{
  NOX::Parameter::List& outputList = params.sublist("Output");

  outputList.setParameter("Nonlinear Iterations", niter);
  outputList.setParameter("2-Norm of Residual", solnptr->getNormF());
  outputList.setParameter("Total Number of Line Search Iterations", 
			  totalNumLSIterations);
  outputList.setParameter("Total Number of Failed Line Searches", 
			  totalNumFailedLineSearches);
  outputList.setParameter("Total Number Steps Requiring Line Search", 
			  totalNumLSSteps);

  return true;
}

bool TensorBased::performStandardLinesearch()
{
  // Copy pointers into temporary references
  Abstract::Group& newsoln = *solnptr;

  // Local variables
  int lsIterations = 1;
  bool ok;
  bool isFailed = false;
  double newf;
  double oldf = 0.5*oldsoln.getNormF()*oldsoln.getNormF();  
                            // Redefined f() = 0.5 * F'*F = 0.5*||F||^2


  // Computation of directional derivative used in curvature condition
  // Note that for Newton direction, fprime = -2.0*oldf
  Abstract::Vector* tmpvecptr = oldsoln.getX().clone(ShapeCopy);
  oldsoln.applyJacobian(dir,*tmpvecptr);
  double fprime = tmpvecptr->dot(oldsoln.getF());

  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Std. Tensor Line Search -- \n";
  }

  // Compute first trial point and its function value
  step = defaultStep;
  newsoln.computeX(oldsoln, dir, step);
  newsoln.computeF();    
  newf = 0.5*newsoln.getNormF()*newsoln.getNormF();  

  // Compute the convergence criteria for the line search 
  double convergence = oldf + alpha*step*fprime;

  // Is a linesearch needed?
  if (newf >= convergence) {

    totalNumLSSteps++;

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << lsIterations << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << endl;
    }
      
    // Is the tensor step a descent direction?  If not, switch to Newton
    if (fprime >= 0) {
      dir = oldsoln.getNewton();
      oldsoln.applyJacobian(dir,*tmpvecptr);
      fprime = tmpvecptr->dot(oldsoln.getF());
    }
    else {
      step = step * 0.5;
    }

    // Update the number of linesearch iterations
    lsIterations ++;
      
    // Compute trial point and its function value
    newsoln.computeX(oldsoln, dir, step);
    newsoln.computeF();
    newf = 0.5*newsoln.getNormF()*newsoln.getNormF(); 

    // Recompute convergence criteria based on new step or new fprime
    convergence = oldf + alpha*step*fprime;

    //Iterate until Armijo-Goldstein condition is satisfied
    while (newf >= convergence) {  
      
      if (Utils::doPrint(Utils::InnerIteration)) {
	cout << setw(3) << lsIterations << ":";
	cout << " step = " << Utils::sci(step);
	cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
	cout << " newf = " << Utils::sci(sqrt(2.*newf));
	cout << endl;
      }
	
      lsIterations ++;
      step = step * 0.5;

      // Check for linesearch failure: if true, exit the method
      if ((step < minStep) || (lsIterations > maxIters)) {

	totalNumFailedLineSearches += 1;
	step = recoveryStep;

	// Update the group and exit linesearch
	newsoln.computeX(oldsoln, dir, step);
	newsoln.computeF();    
	newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

	// Print out info
	cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
	cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
	cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
	cout << " (USING RECOVERY STEP!)" << endl;
	cout << Utils::fill(72) << "\n" << endl;

	isFailed = true;
	return (!isFailed);
      }

      // Compute new trial point and its function value
      newsoln.computeX(oldsoln, dir, step);
      newsoln.computeF();    
      newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

      // Recompute convergence criteria based on new step
      convergence = oldf + alpha*step*fprime;

    }
  }

  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << setw(3) << lsIterations << ":";
    cout << " step = " << Utils::sci(step);
    cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
    cout << " newf = " << Utils::sci(sqrt(2.*newf));
    if (isFailed)
      cout << " (USING RECOVERY STEP!)" << endl;
    else
      cout << " (STEP ACCEPTED!)" << endl;
    cout << Utils::fill(72) << "\n" << endl;
  }
  
  ok = !isFailed;
  totalNumLSIterations += lsIterations - 1;

  delete tmpvecptr;

  return ok;
}


bool TensorBased::performCurvilinearLinesearch()
{
  // Copy pointer into temporary reference
  Abstract::Group& newsoln = *solnptr;

  // local variables
  int lsIterations = 1;
  bool ok;
  bool isFailed = false;
  double newf;
  double oldf = 0.5*oldsoln.getNormF()*oldsoln.getNormF();  
                            // Redefined f() = 0.5 * F'*F = 0.5*||F||^2

  // Computation of directional derivative used in curvature condition
  Abstract::Vector* tmpvecptr = oldsoln.getX().clone(ShapeCopy);
  oldsoln.applyJacobian(oldsoln.getNewton(),*tmpvecptr);   // bwb - dNewton 
  double fprime = tmpvecptr->dot(oldsoln.getF());
  
  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Curvilinear Linesearch -- \n";
  }
  
  // Compute curvilinear step, new trial point, and its function value
  step = defaultStep;
  newsoln.computeX(oldsoln, dir, step);
  newsoln.computeF();    
  newf = 0.5*newsoln.getNormF()*newsoln.getNormF();  
  
  // Compute the convergence criteria for the line search 
  double convergence = oldf + alpha*step*fprime;

  // Increment the counter of Newton steps requiring a linesearch
  if (newf >= convergence)
    totalNumLSSteps++;

  //Iterate until Armijo-Goldstein condition is satisfied
  while (newf >= convergence) {  
    
    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << lsIterations << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << endl;
    }
    
    lsIterations ++;

    // Backtrack by halving the step
    step = step * 0.5;

    // Check for linesearch failure: if true, exit the method
    if ((step < minStep) || (lsIterations > maxIters)) {

      totalNumFailedLineSearches += 1;
      step = recoveryStep;

      // Update the group, and exit linesearch
      ok = direction.computeCurvilinearStep(dir, newsoln, *this, step);
      newsoln.computeX(oldsoln, dir, 1.0);
      newsoln.computeF();    
      newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

      // Print out info
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (USING RECOVERY STEP!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;

      isFailed = true;
      return (!isFailed);
    }

    // Compute curvilinear step, new trial point, and its function value
    ok = direction.computeCurvilinearStep(dir, newsoln, *this, step);   
    newsoln.computeX(oldsoln, dir, 1.0);
    newsoln.computeF();    
    newf = 0.5*newsoln.getNormF()*newsoln.getNormF();

    // Recompute convergence criteria based on new step
    convergence = oldf + alpha*step*fprime;

  }  // end while loop


  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << setw(3) << lsIterations << ":";
    cout << " step = " << Utils::sci(step);
    cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
    cout << " newf = " << Utils::sci(sqrt(2.*newf));
    cout << " (STEP ACCEPTED!)" << endl;
    cout << Utils::fill(72) << "\n" << endl;
  }

#if DEBUG_LEVEL > 1
  cout << "Current point: ";
  newsoln.getX().print();
#endif

  ok = !isFailed;
  totalNumLSIterations += lsIterations - 1;

  delete tmpvecptr;

  return ok;
}


#endif
