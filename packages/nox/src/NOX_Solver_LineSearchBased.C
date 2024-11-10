// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Solver_LineSearchBased.H"    // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Observer.hpp"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_LineSearch_Generic.H"
#include "NOX_LineSearch_Factory.H"
#include "NOX_Direction_Generic.H"
#include "NOX_Direction_Factory.H"
#include "NOX_SolverStats.hpp"

NOX::Solver::LineSearchBased::
LineSearchBased(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
        const Teuchos::RCP<NOX::StatusTest::Generic>& t,
        const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(xGrp),                               // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  dirPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  testPtr(t),
  paramsPtr(p)
{
  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils(); 
  observer = NOX::Solver::parseObserver(p->sublist("Solver Options"));
  this->init();
}

// Protected
void NOX::Solver::LineSearchBased::init()
{
  // Initialize
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  lineSearchPtr = NOX::LineSearch::
    buildLineSearch(globalDataPtr, paramsPtr->sublist("Line Search"));

  directionPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Direction"));

  if ( paramsPtr->isType<bool>( "Catch Throws During Solve" ) )
    catchThrowsDuringSolve = paramsPtr->get<bool>( "Catch Throws During Solve" );
  else
    catchThrowsDuringSolve = false;

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

}

void NOX::Solver::LineSearchBased::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initialGuess);
  testPtr = t;
  globalDataPtr->getNonConstSolverStatistics()->reset();
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
}

void NOX::Solver::LineSearchBased::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  globalDataPtr->getNonConstSolverStatistics()->reset();
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
}

void NOX::Solver::LineSearchBased::
reset()
{
  globalDataPtr->getNonConstSolverStatistics()->reset();
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
}

NOX::Solver::LineSearchBased::~LineSearchBased()
{

}


NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::getStatus() const
{
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::step()
{
  observer->runPreIterate(*this);

  // On the first step, do some initializations
  if (nIter == 0) {
    globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearSolves();

    // Compute F of initital guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      utilsPtr->out() << "NOX::Solver::LineSearchBased::init - "
              << "Unable to compute F" << std::endl;
      throw std::runtime_error("NOX Error");
    }

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) &&
    (utilsPtr->isPrintType(NOX::Utils::Warning))) {
      utilsPtr->out() << "Warning: NOX::Solver::LineSearchBased::init() - "
              << "The solution passed into the solver (either "
              << "through constructor or reset method) "
              << "is already converged!  The solver wil not "
              << "attempt to solve this system since status is "
              << "flagged as converged." << std::endl;
    }

    if (status == NOX::StatusTest::Unconverged) printUpdate();
  }

  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = directionPtr->compute(*dirPtr, soln, *this);
  if (!ok)
  {
    utilsPtr->out() << "NOX::Solver::LineSearchBased::iterate - unable to calculate direction" << std::endl;
    status = NOX::StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Update iteration count.
  nIter ++;
  globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearIterations();

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;

  // Do line search and compute new soln.
  observer->runPreSolutionUpdate(*dirPtr,*this);
  observer->runPreLineSearch(*this);
  ok = lineSearchPtr->compute(soln, stepSize, *dirPtr, *this);
  observer->runPostLineSearch(*this);
  observer->runPostSolutionUpdate(*this);
  if (!ok)
  {
    if (stepSize == 0.0)
    {
      utilsPtr->out() << "NOX::Solver::LineSearchBased::iterate - line search failed" << std::endl;
      status = NOX::StatusTest::Failed;
      observer->runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out() << "NOX::Solver::LineSearchBased::iterate - using recovery step for line search" << std::endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utilsPtr->out() << "NOX::Solver::LineSearchBased::iterate - unable to compute F" << std::endl;
    status = NOX::StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Evaluate the current status.
  status = test.checkStatus(*this, checkType);

  observer->runPostIterate(*this);

  printUpdate();

  return status;
}

NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::solve()
{
  observer->runPreSolve(*this);

  this->reset();

  // Iterate until converged or failed
  if ( catchThrowsDuringSolve )
  {
    try
    {
      while (status == NOX::StatusTest::Unconverged)
        step();
    }
    catch ( std::exception & e )
    {
      utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
      utilsPtr->out()
        << "-- WARNING --\n"
        << "NOX::Solver::LineSearchBased::solve caught exception during solver iteration:\n\n"
        << e.what() << "\n\n"
        << "Setting solver status to Failed and continuing anyway since catchThrowsDuringSolve is true.\n";
      utilsPtr->out() << NOX::Utils::fill(72) << "\n" << std::endl;

      status = NOX::StatusTest::Failed;
    }
  }
  else
  {
    while (status == NOX::StatusTest::Unconverged)
      step();
  }

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  observer->runPostSolve(*this);

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
  return *oldSolnPtr;
}

int NOX::Solver::LineSearchBased::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList&
NOX::Solver::LineSearchBased::getList() const
{
  return *paramsPtr;
}

Teuchos::RCP<const NOX::SolverStats>
NOX::Solver::LineSearchBased::getSolverStatistics() const
{
  return globalDataPtr->getSolverStatistics();
}

// protected
void NOX::Solver::LineSearchBased::printUpdate()
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested
  if ((status == NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dirPtr->norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "||F|| = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}

double NOX::Solver::LineSearchBased::getStepSize() const
{
  return stepSize;
}
