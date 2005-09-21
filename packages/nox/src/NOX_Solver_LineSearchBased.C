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
#include "NOX_Utils.H"

NOX::Solver::LineSearchBased::
LineSearchBased(const Teuchos::RefCountPtr<NOX::Abstract::Group>& xGrp, 
		const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
		const Teuchos::RefCountPtr<NOX::Parameter::List>& p) :
  solnPtr(xGrp),		        // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  oldSoln(*oldSolnPtr),		        // reference to just-created pointer
  dirPtr(xGrp->getX().clone(ShapeCopy)), // create via clone 
  dir(*dirPtr),			        // reference to just-created pointer
  testPtr(t),			// pointer to t
  paramsPtr(p),		// pointer to p
  utils(paramsPtr->sublist("Printing")),                // intialize the utils
  lineSearch(utils, paramsPtr->sublist("Line Search")), // initialize the line search
  direction(utils, paramsPtr->sublist("Direction")),     // initialize the direction
  prePostOperator(utils, paramsPtr->sublist("Solver Options"))
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

  // Get the checktype
  checkType = (NOX::StatusTest::CheckType) paramsPtr->
    sublist("Solver Options").getParameter("Status Test Check Type", 
					   NOX::StatusTest::Minimal);

  // Print out parameters
  if (utils.isPrintType(NOX::Utils::Parameters)) 
  {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils.out(),5);
  }

  // Compute F of initital guess
  NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    utils.out() << "NOX::Solver::LineSearchBased::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testPtr->checkStatus(*this, checkType);
  if ((status == NOX::StatusTest::Converged) &&
      (utils.isPrintType(NOX::Utils::Warning)))
  {
    utils.out() << "Warning: NOX::Solver::LineSearchBased::init() - The solution passed "
	 << "into the solver (either through constructor or reset method) "
	 << "is already converged!  The solver wil not "
	   << "attempt to solve this system since status is flagged as "
	 << "converged." << endl;
  }

  // Print out status tests
  if (utils.isPrintType(NOX::Utils::Parameters)) 
  {
    utils.out() << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(utils.out(), 5);
    utils.out() <<"\n" << NOX::Utils::fill(72) << "\n";
  }

}

bool NOX::Solver::LineSearchBased::
reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& xGrp, 
      const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
      const Teuchos::RefCountPtr<NOX::Parameter::List>& p) 
{
  solnPtr = xGrp;
  testPtr = t;
  paramsPtr = p;		
  utils.reset(paramsPtr->sublist("Printing"));
  lineSearch.reset(paramsPtr->sublist("Line Search"));	
  direction.reset(paramsPtr->sublist("Direction"));
  prePostOperator.reset(utils, paramsPtr->sublist("Solver Options"));

  init();

  return true;
}

bool NOX::Solver::LineSearchBased::
reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& xGrp, 
      const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t)
{
  solnPtr = xGrp;
  testPtr = t;
  init();
  return true;
}

NOX::Solver::LineSearchBased::~LineSearchBased() 
{
  delete oldSolnPtr;
  delete dirPtr;
}


NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::iterate()
{
  prePostOperator.runPreIterate(*this);

  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    prePostOperator.runPostIterate(*this);
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
    utils.out() << "NOX::Solver::LineSearchBased::iterate - unable to calculate direction" << endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
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
      utils.out() << "NOX::Solver::LineSearchBased::iterate - line search failed" << endl;
      status = NOX::StatusTest::Failed;
      prePostOperator.runPostIterate(*this);
      return status;
    }
    else if (utils.isPrintType(NOX::Utils::Warning))
      utils.out() << "NOX::Solver::LineSearchBased::iterate - using recovery step for line search" << endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    utils.out() << "NOX::Solver::LineSearchBased::iterate - unable to compute F" << endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }

  // Evaluate the current status.
  status = test.checkStatus(*this, checkType);
 
  prePostOperator.runPostIterate(*this);

  // Return status.
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::solve()
{
  prePostOperator.runPreSolve(*this);

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

  prePostOperator.runPostSolve(*this);

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
      (utils.isPrintType(NOX::Utils::OuterIterationStatusTest))) 
  {
    utils.out() << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Status Test Results --\n";    
    testPtr->print(utils.out());
    utils.out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utils.isPrintType(NOX::Utils::OuterIteration)) 
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utils.isPrintType(NOX::Utils::OuterIteration)) 
  {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utils.out() << "f = " << utils.sciformat(normSoln);
    utils.out() << "  step = " << utils.sciformat(step);
    utils.out() << "  dx = " << utils.sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      utils.out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utils.out() << " (Failed!)";
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) && 
      (utils.isPrintType(NOX::Utils::OuterIteration))) 
  {
    utils.out() << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Final Status Test Results --\n";    
    testPtr->print(utils.out());
    utils.out() << NOX::Utils::fill(72) << "\n";
  }
}

double NOX::Solver::LineSearchBased::getStepSize() const
{
  return step;
}
