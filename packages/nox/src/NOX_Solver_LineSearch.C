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

#include "NOX_Solver_LineSearch.H"	// class definition
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

LineSearch::LineSearch(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p) :
  solnptr(&xgrp),		// pointer to xgrp
  oldsolnptr(xgrp.clone(DeepCopy)), // create via clone
  oldsoln(*oldsolnptr),		// reference to just-created pointer
  dirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  dir(*dirptr),			// reference to just-created pointer
  testptr(&t),			// pointer to t
  iparams(p),			// copy p
  linesearch(iparams.sublist("Line Search")), // initialize line search
  direction(iparams.sublist("Direction")) // initialize direction
{
  init();
}

// Protected
void LineSearch::init()
{
  // Initialize 
  step = 0;
  niter = 0;
  status = Status::Unconverged;

  // Set up utilities (i.e., set print processor, etc)
  Utils::setUtils(iparams);
  
  // Print out initialization information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    iparams.print(cout,5);
  }

  // Compute RHS of initital guess
  bool ok = solnptr->computeRHS();
  if (!ok) {
    cout << "NOX::Solver::LineSearch::init - Unable to compute RHS" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testptr->operator()(*this);
  if (status == Status::Converged) {
    if (Utils::doPrint(Utils::Warning)) {
      cout << "Warning: NOX::Solver::LineSearch::init() - The solution passed "
	   << "into the solver (either through constructor or reset method) "
	   << "is already converged!  The solver wil not "
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

bool LineSearch::reset(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p) 
{
  solnptr = &xgrp;
  testptr = &t;
  iparams = p;			
  linesearch.reset(iparams.sublist("Line Search"));	
  direction.reset(iparams.sublist("Direction"));
  init();
  return true;
}

LineSearch::~LineSearch() 
{
  delete oldsolnptr;
  delete dirptr;
}


Status::StatusType LineSearch::getStatus()
{
  return status;
}

Status::StatusType LineSearch::iterate()
{
  // First check status
  if (status != Status::Unconverged) 
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnptr;
  Status::Test& test = *testptr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction(dir, soln, *this);
  if (!ok) {
    cout << "NOX::Solver::LineSearch::iterate - unable to calculate direction" << endl;
    status = Status::Failed;
    return status;
  }

  // Copy current soln to the old soln.
  oldsoln = soln;

  // Do line search and compute new soln.
  ok = linesearch(soln, step, oldsoln, dir);
  if (!ok) {
    if (step == 0) {
      cout << "NOX::Solver::LineSearch::iterate - linesearch failed" << endl;
      status = Status::Failed;
      return status;
    }
    else if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Solver::LineSearch::iterate - using recovery step for linesearch" << endl;
  }

  // Compute RHS for new current solution.
  ok = soln.computeRHS();
  if (!ok) {
    cout << "NOX::Solver::LineSearch::iterate - unable to compute RHS" << endl;
    status = Status::Failed;
    return status;
  }

  // Update iteration count.
  niter ++;

  // Evaluate the current status.
  status = test(*this);
 
  // Return status.
  return status;
}

Status::StatusType LineSearch::solve()
{
  printUpdate();

  // Iterate until converged or failed
  while (status == Status::Unconverged) {
    status = iterate();
    printUpdate();
  }

  return status;
}

const Abstract::Group& LineSearch::getSolutionGroup() const
{
  return *solnptr;
}

const Abstract::Group& LineSearch::getPreviousSolutionGroup() const
{
  return oldsoln;
}

int LineSearch::getNumIterations() const
{
  return niter;
}

const Parameter::List& LineSearch::getOutputParameters() const
{
  oparams.setParameter("Nonlinear Iterations", niter);
  oparams.setParameter("2-Norm of Residual", solnptr->getNormRHS());
  return oparams;
}

// protected
void LineSearch::printUpdate() 
{
  double norm_soln;
  double norm_step;

  // All processes participate in the computation of these norms...
  if (Utils::doAllPrint(Utils::OuterIteration)) {
    norm_soln = solnptr->getNormRHS();
    norm_step = (niter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (Utils::doPrint(Utils::OuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Line Search Method Step " << niter << " -- \n";
    cout << "f = " << Utils::sci(norm_soln);
    cout << "  step = " << Utils::sci(step);
    cout << "  dx = " << Utils::sci(norm_step);
    if (status == Status::Converged)
      cout << " (Converged!)";
    if (status == Status::Failed)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != Status::Unconverged) && 
      (Utils::doPrint(Utils::OuterIteration))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testptr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}

