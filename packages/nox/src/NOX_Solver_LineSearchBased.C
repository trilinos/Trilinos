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

#include "NOX_Solver_LineSearchBased.H"	// class definition
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

LineSearchBased::LineSearchBased(Abstract::Group& xGrp, 
				 StatusTest::Generic& t, 
				 const Parameter::List& p) :
  solnPtr(&xGrp),		// pointer to xGrp
  oldSolnPtr(xGrp.clone(DeepCopy)), // create via clone
  oldSoln(*oldSolnPtr),		// reference to just-created pointer
  dirPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  dir(*dirPtr),			// reference to just-created pointer
  testPtr(&t),			// pointer to t
  params(p),			// copy p
  lineSearch(params.sublist("Line Search")), // initialize the line search
  direction(params.sublist("Direction")) // initialize the direction
{
  init();
}

// Protected
void LineSearchBased::init()
{
  // Initialize 
  step = 0;
  nIter = 0;
  status = StatusTest::Unconverged;

  // Set up utilities (i.e., set print processor, etc)
  Utils::setUtils(params);
  
  // Print out initialization information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    params.print(cout,5);
  }

  // Compute F of initital guess
  bool ok = solnPtr->computeF();
  if (!ok) {
    cout << "NOX::Solver::LineSearchBased::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testPtr->checkStatus(*this);
  if (status == StatusTest::Converged) {
    if (Utils::doPrint(Utils::Warning)) {
      cout << "Warning: NOX::Solver::LineSearchBased::init() - The solution passed "
	   << "into the solver (either through constructor or reset method) "
	   << "is already converged!  The solver wil not "
	   << "attempt to solve this system since status is flagged as "
	   << "converged." << endl;
    }
  }

  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";
  }

}

bool LineSearchBased::reset(Abstract::Group& xGrp, StatusTest::Generic& t, const Parameter::List& p) 
{
  solnPtr = &xGrp;
  testPtr = &t;
  params = p;			
  lineSearch.reset(params.sublist("Line Search"));	
  direction.reset(params.sublist("Direction"));
  init();
  return true;
}

LineSearchBased::~LineSearchBased() 
{
  delete oldSolnPtr;
  delete dirPtr;
}


NOX::StatusTest::StatusType LineSearchBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType LineSearchBased::iterate()
{
  // First check status
  if (status != StatusTest::Unconverged) 
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnPtr;
  StatusTest::Generic& test = *testPtr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction.compute(dir, soln, *this);
  if (!ok) {
    cout << "NOX::Solver::LineSearchBased::iterate - unable to calculate direction" << endl;
    status = StatusTest::Failed;
    return status;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  oldSoln = soln;

  // Do line search and compute new soln.
  ok = lineSearch.compute(soln, step, dir, *this);
  if (!ok) {
    if (step == 0) {
      cout << "NOX::Solver::LineSearchBased::iterate - line search failed" << endl;
      status = StatusTest::Failed;
      return status;
    }
    else if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Solver::LineSearchBased::iterate - using recovery step for line search" << endl;
  }

  // Compute F for new current solution.
  ok = soln.computeF();
  if (!ok) {
    cout << "NOX::Solver::LineSearchBased::iterate - unable to compute F" << endl;
    status = StatusTest::Failed;
    return status;
  }

  // Evaluate the current status.
  status = test.checkStatus(*this);
 
  // Return status.
  return status;
}

NOX::StatusTest::StatusType LineSearchBased::solve()
{
  printUpdate();

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = iterate();
    printUpdate();
  }

  Parameter::List& outputParams = params.sublist("Output");
  outputParams.setParameter("Nonlinear Iterations", nIter);
  outputParams.setParameter("2-Norm of Residual", solnPtr->getNormF());

  return status;
}

const Abstract::Group& LineSearchBased::getSolutionGroup() const
{
  return *solnPtr;
}

const Abstract::Group& LineSearchBased::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int LineSearchBased::getNumIterations() const
{
  return nIter;
}

const Parameter::List& LineSearchBased::getParameterList() const
{
  return params;
}

// protected
void LineSearchBased::printUpdate() 
{
  double normSoln;
  double normStep;

  // Print the status test parameters at each iteration if requested  
  if ((status == StatusTest::Unconverged) && 
      (Utils::doPrint(Utils::OuterIterationStatusTest))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Status Test Results --\n";    
    testPtr->print(cout);
    cout << Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (Utils::doAllPrint(Utils::OuterIteration)) {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (Utils::doPrint(Utils::OuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Nonlinear Solver Step " << nIter << " -- \n";
    cout << "f = " << Utils::sci(normSoln);
    cout << "  step = " << Utils::sci(step);
    cout << "  dx = " << Utils::sci(normStep);
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
    testPtr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}


double LineSearchBased::getStepSize() const
{
  return step;
}
