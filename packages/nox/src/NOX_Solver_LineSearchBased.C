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

#include "NOX_Solver_LineSearchBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Parameter_PrePostOperator.H"
#include "NOX_Utils.H"

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

NOX::Solver::LineSearchBased::LineSearchBased(NOX::Abstract::Group& xGrp, 
					      NOX::StatusTest::Generic& t, 
					      NOX::Parameter::List& p) :
  solnPtr(&xGrp),		        // pointer to xGrp
  oldSolnPtr(xGrp.clone(DeepCopy)),     // create via clone
  oldSoln(*oldSolnPtr),		        // reference to just-created pointer
  dirPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  dir(*dirPtr),			        // reference to just-created pointer
  testPtr(&t),			// pointer to t
  paramsPtr(&p),		// pointer to p
  utils(paramsPtr->sublist("Printing")),                // intialize the utils
  lineSearch(utils, paramsPtr->sublist("Line Search")), // initialize the line search
  direction(utils, paramsPtr->sublist("Direction")),     // initialize the direction
  prePostOperatorPtr(0),
  havePrePostOperator(false)
{
  init();
}

// Protected
void NOX::Solver::LineSearchBased::init()
{
  // Initialize 
  step = 0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;

  // Check for a user defined Pre/Post Operator
  NOX::Parameter::List& p = paramsPtr->sublist("Solver Options");
  havePrePostOperator = false;
  prePostOperatorPtr = 0;
  if (p.isParameter("User Defined Pre/Post Operator")) {
    if (p.isParameterArbitrary("User Defined Pre/Post Operator")) {
      prePostOperatorPtr = dynamic_cast<NOX::Parameter::PrePostOperator*>
	(p.getArbitraryParameter("User Defined Pre/Post Operator").clone());
      if (prePostOperatorPtr != 0)
	havePrePostOperator = true;
      else
	if (utils.isPrintProcessAndType(NOX::Utils::Warning))
	  cout << "Warning: NOX::Solver::LineSearchBased::init() - " 
	       << "\"User Defined Pre/Post Operator\" not derived from " 
	       << "NOX::Parameter::PrePostOperator class!\n" 
	       << "Ignoring this flag!"<< endl;
    }
    else {
      cout << "ERROR: NOX::Solver::LineSearchBased::init() - the parameter "
	   << "\"User Defined Pre/Post Operator\" must be derived from an"
	   << "arbitrary parameter!" << endl;
      throw "NOX Error";
    }
  }
  
  // Print out parameters
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) 
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(cout,5);
  }

  // Compute F of initital guess
  NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    cout << "NOX::Solver::LineSearchBased::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testPtr->checkStatus(*this);
  if ((status == NOX::StatusTest::Converged) &&
      (utils.isPrintProcessAndType(NOX::Utils::Warning)))
  {
    cout << "Warning: NOX::Solver::LineSearchBased::init() - The solution passed "
	 << "into the solver (either through constructor or reset method) "
	 << "is already converged!  The solver wil not "
	   << "attempt to solve this system since status is flagged as "
	 << "converged." << endl;
  }

  // Print out status tests
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) 
  {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << NOX::Utils::fill(72) << "\n";
  }

}

bool NOX::Solver::LineSearchBased::reset(NOX::Abstract::Group& xGrp, 
					 NOX::StatusTest::Generic& t, 
					 NOX::Parameter::List& p) 
{
  solnPtr = &xGrp;
  testPtr = &t;
  paramsPtr = &p;		
  utils.reset(paramsPtr->sublist("Printing"));
  lineSearch.reset(paramsPtr->sublist("Line Search"));	
  direction.reset(paramsPtr->sublist("Direction"));
  init();
  return true;
}

bool NOX::Solver::LineSearchBased::reset(NOX::Abstract::Group& xGrp, 
					 NOX::StatusTest::Generic& t)
{
  solnPtr = &xGrp;
  testPtr = &t;
  init();
  return true;
}

NOX::Solver::LineSearchBased::~LineSearchBased() 
{
  delete prePostOperatorPtr;
  delete oldSolnPtr;
  delete dirPtr;
}


NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::iterate()
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreIterate(*this);

  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction.compute(dir, soln, *this);
  if (!ok) 
  {
    cout << "NOX::Solver::LineSearchBased::iterate - unable to calculate direction" << endl;
    status = NOX::StatusTest::Failed;
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  oldSoln = soln;

  // Do line search and compute new soln.
  ok = lineSearch.compute(soln, step, dir, *this);
  if (!ok) 
  {
    if (step == 0) 
    {
      cout << "NOX::Solver::LineSearchBased::iterate - line search failed" << endl;
      status = NOX::StatusTest::Failed;
      if (havePrePostOperator)
	prePostOperatorPtr->runPostIterate(*this);
      return status;
    }
    else if (utils.isPrintProcessAndType(NOX::Utils::Warning))
      cout << "NOX::Solver::LineSearchBased::iterate - using recovery step for line search" << endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    cout << "NOX::Solver::LineSearchBased::iterate - unable to compute F" << endl;
    status = NOX::StatusTest::Failed;
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  // Evaluate the current status.
  status = test.checkStatus(*this);
 
  if (havePrePostOperator)
    prePostOperatorPtr->runPostIterate(*this);

  // Return status.
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::solve()
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreSolve(*this);

  printUpdate();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged) 
  {
    status = iterate();
    printUpdate();
  }

  NOX::Parameter::List& outputParams = paramsPtr->sublist("Output");
  outputParams.setParameter("Nonlinear Iterations", nIter);
  outputParams.setParameter("2-Norm of Residual", solnPtr->getNormF());

  if (havePrePostOperator)
    prePostOperatorPtr->runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group& 
NOX::Solver::LineSearchBased::getSolutionGroup() const
{
  return *solnPtr;
}

const NOX::Abstract::Group& 
NOX::Solver::LineSearchBased::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int NOX::Solver::LineSearchBased::getNumIterations() const
{
  return nIter;
}

const NOX::Parameter::List& 
NOX::Solver::LineSearchBased::getParameterList() const
{
  return *paramsPtr;
}

// protected
void NOX::Solver::LineSearchBased::printUpdate() 
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested  
  if ((status == NOX::StatusTest::Unconverged) && 
      (utils.isPrintProcessAndType(NOX::Utils::OuterIterationStatusTest))) 
  {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Status Test Results --\n";    
    testPtr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utils.isPrintType(NOX::Utils::OuterIteration)) 
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utils.isPrintProcessAndType(NOX::Utils::OuterIteration)) 
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Nonlinear Solver Step " << nIter << " -- \n";
    cout << "f = " << utils.sciformat(normSoln);
    cout << "  step = " << utils.sciformat(step);
    cout << "  dx = " << utils.sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      cout << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      cout << " (Failed!)";
    cout << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) && 
      (utils.isPrintProcessAndType(NOX::Utils::OuterIteration))) 
  {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testPtr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }
}

double NOX::Solver::LineSearchBased::getStepSize() const
{
  return step;
}
