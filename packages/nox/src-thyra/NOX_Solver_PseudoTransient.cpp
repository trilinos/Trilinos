// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Solver_PseudoTransient.hpp"    // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_LineSearch_Generic.H"
#include "NOX_LineSearch_Factory.H"
#include "NOX_Direction_Generic.H"
#include "NOX_Direction_Factory.H"
#include "NOX_Observer.hpp"
#include "NOX_SolverStats.hpp"
#include <limits>

// This relies explicitly on thyra
#include "NOX_Thyra_Group.H"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ModelEvaluatorBase.hpp"

NOX::Solver::PseudoTransient::
PseudoTransient(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
        const Teuchos::RCP<NOX::StatusTest::Generic>& t,
        const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(xGrp),                               // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  transientResidualGroup(xGrp->clone(DeepCopy)),     // create via clone
  dirPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  testPtr(t)
{
  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils();
  observer = NOX::Solver::parseObserver(p->sublist("Solver Options"));

  this->setMyParamList(p);
  thyraSolnGroup = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>(solnPtr,true);
  thyraOldSolnGroup = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>(oldSolnPtr,true);
  thyraTransientResidualGroup = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>(transientResidualGroup,true);

  init();
}

// Protected
void NOX::Solver::PseudoTransient::init()
{
  // Initialize
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();

  Teuchos::RCP<Teuchos::ParameterList> paramsPtr = this->getMyNonconstParamList();
  paramsPtr->validateParametersAndSetDefaults(*this->getValidParameters());

  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  lineSearchPtr = NOX::LineSearch::
    buildLineSearch(globalDataPtr, paramsPtr->sublist("Line Search"));

  directionPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Direction"));

  deltaInit = paramsPtr->sublist("Pseudo-Transient").get<double>("deltaInit");
  delta = deltaInit;
  inv_delta = 1.0 / delta;
  deltaMax = paramsPtr->sublist("Pseudo-Transient").get<double>("deltaMax");
  deltaMin = paramsPtr->sublist("Pseudo-Transient").get<double>("deltaMin");
  time = 0.0;

  use_transient_residual =
    paramsPtr->sublist("Pseudo-Transient").get<bool>("Use Transient Residual in Direction Computation");

  max_pseudo_transient_iterations =
    paramsPtr->sublist("Pseudo-Transient").get<int>("Maximum Number of Pseudo-Transient Iterations");

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

}

void NOX::Solver::PseudoTransient::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
{
  this->setMyParamList(paramList);
}

Teuchos::RCP<const Teuchos::ParameterList> NOX::Solver::PseudoTransient::getValidParameters() const
{
  if (is_null(validParameters)) {
    validParameters = Teuchos::parameterList();

    validParameters->sublist("Solver Options").disableRecursiveValidation();
    validParameters->sublist("Line Search").disableRecursiveValidation();
    validParameters->sublist("Direction").disableRecursiveValidation();
    validParameters->sublist("Printing").disableRecursiveValidation();
    validParameters->sublist("Output").disableRecursiveValidation();
    validParameters->sublist("Thyra Group Options").disableRecursiveValidation();
    validParameters->sublist("Status Tests").disableRecursiveValidation();

    validParameters->set<std::string>("Nonlinear Solver","Pseudo-Transient");

    Teuchos::ParameterList& ps_params = validParameters->sublist("Pseudo-Transient");

    ps_params.set<double>("deltaInit",1.0e-4,"Initial time step size.");
    ps_params.set<double>("deltaMax",1.0e+0,"Maximum time step size.  If the new step size is greater than this value, the transient terms will be eliminated from the Newton interation resulting in a full Newton solve.");
    ps_params.set<double>("deltaMin",1.0e-5, "Minimum step size.");
    ps_params.set<bool>("Use Transient Residual in Direction Computation",true);
    ps_params.set<int>("Maximum Number of Pseudo-Transient Iterations", std::numeric_limits<int>::max());

  }

  return validParameters;
}

void NOX::Solver::PseudoTransient::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initialGuess);
  testPtr = t;
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();
  time = 0.0;
}

void NOX::Solver::PseudoTransient::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();
  time = 0.0;
}

void NOX::Solver::PseudoTransient::
reset()
{
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();
  time = 0.0;
}

NOX::StatusTest::StatusType
NOX::Solver::PseudoTransient::getStatus() const
{
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::PseudoTransient::step()
{
  observer->runPreIterate(*this);

  // On the first step, do some initializations
  if (nIter == 0) {
    globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearSolves();

    // Stupid row sum scaling.  Computes a Jacobian before x_dot or
    // alpha/beta is known so first Jacobian is for steady-state.  Need
    // to reset this as we enter the solver.
    solnPtr->setX(solnPtr->getX());

    // Compute F of initital guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      utilsPtr->out() << "NOX::Solver::PseudoTransient::init - "
              << "Unable to compute F" << std::endl;
      throw std::runtime_error("NOX Error");
    }

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) &&
    (utilsPtr->isPrintType(NOX::Utils::Warning))) {
      utilsPtr->out() << "Warning: NOX::Solver::PseudoTransient::init() - "
              << "The solution passed into the solver (either "
              << "through constructor or reset method) "
              << "is already converged!  The solver wil not "
              << "attempt to solve this system since status is "
              << "flagged as converged." << std::endl;
    }

    printUpdate();
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

  // Pseudo-transient: change the Jacobian evaluation to evaluate a transient version
  if (nIter < max_pseudo_transient_iterations) {
    deltaOld = delta;

    // Update step size
    if (nIter == 0)
      delta = deltaInit;
    else
      delta = deltaOld * oldSolnPtr->getNormF() / solnPtr->getNormF();

    inv_delta = 1.0 / delta;
    if (delta > deltaMax)
      inv_delta = 0.0;

    if (delta < deltaMin) {
      delta = deltaMin;
      inv_delta = 1.0 / delta;
    }

    time += delta;

    Teuchos::RCP<const ::Thyra::VectorBase<double> > x =
      thyraSolnGroup->get_current_x();

    Teuchos::RCP<const ::Thyra::VectorBase<double> > x_old =
      thyraOldSolnGroup->get_current_x();

    // Compute x_dot using forward finite difference
    if (is_null(x_dot))
      x_dot = ::Thyra::createMember(x->space());

    if (nIter == 0)
      ::Thyra::put_scalar(0.0,x_dot.ptr());
    else {
      ::Thyra::V_StVpStV(x_dot.ptr(),inv_delta,*x,-inv_delta,*x_old);
    }

    thyraSolnGroup->enablePseudoTransientTerms(x_dot,inv_delta,1.0,time);
  }
  else {
    delta = std::numeric_limits<double>::max();
    inv_delta = 0.0;
  }

  // Compute the direction for the update vector at the current
  // solution.  Steady-state F is already computed so the only thing
  // to compute is J using our augmented inargs.  If
  // use_transient_residual is true, then we will also compute this
  // quantity and use it.
  bool ok = true;
  if ( use_transient_residual && (nIter < max_pseudo_transient_iterations) ) {
    thyraTransientResidualGroup->setX(solnPtr->getX());
    thyraTransientResidualGroup->enablePseudoTransientTerms(x_dot,inv_delta,1.0,time);
    ok = directionPtr->compute(*dirPtr, *thyraTransientResidualGroup, *this);
    NOX::Abstract::Group::ReturnType rtype = thyraTransientResidualGroup->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      utilsPtr->out() << "NOX::Solver::PseudoTransient::init - "
              << "Unable to compute F" << std::endl;
      throw std::runtime_error("NOX Error");
    }
    thyraTransientResidualGroup->disablePseudoTransientTerms();
  }
  else {
    ok = directionPtr->compute(*dirPtr, soln, *this);
  }

  if (!ok)
  {
    utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - unable to calculate direction" << std::endl;
    status = NOX::StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // reset the inargs to the correct value for steady state residual evaluations
  thyraSolnGroup->disablePseudoTransientTerms();

  // Update iteration count.
  nIter ++;
  globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearIterations();

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;

  // Do line search and compute new soln.
  observer->runPreSolutionUpdate(*dirPtr,*this);
  ok = lineSearchPtr->compute(soln, stepSize, *dirPtr, *this);
  observer->runPostSolutionUpdate(*this);  
  if (!ok)
  {
    if (stepSize == 0.0)
    {
      utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - line search failed" << std::endl;
      status = NOX::StatusTest::Failed;
      observer->runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - using recovery step for line search"
              << std::endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - unable to compute F" << std::endl;
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

NOX::StatusTest::StatusType NOX::Solver::PseudoTransient::solve()
{
  observer->runPreSolve(*this);

  this->reset();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged)
    step();

  Teuchos::ParameterList& outputParams = this->getMyNonconstParamList()->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  observer->runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group&
NOX::Solver::PseudoTransient::getSolutionGroup() const
{
  return *solnPtr;
}

const NOX::Abstract::Group&
NOX::Solver::PseudoTransient::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

int NOX::Solver::PseudoTransient::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList&
NOX::Solver::PseudoTransient::getList() const
{
  return *this->getMyParamList();
}

// protected
void NOX::Solver::PseudoTransient::printUpdate()
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
    utilsPtr->out() << "  delta = " << utilsPtr->sciformat(delta);
    utilsPtr->out() << "  inv_delta = " << utilsPtr->sciformat(inv_delta);
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

double NOX::Solver::PseudoTransient::getStepSize() const
{
  return stepSize;
}

Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::PseudoTransient::getSolutionGroupPtr() const
{return solnPtr;}

Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::PseudoTransient::getPreviousSolutionGroupPtr() const
{return oldSolnPtr;}

Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::PseudoTransient::getListPtr() const
{return this->getMyParamList();}

Teuchos::RCP<const NOX::SolverStats>
NOX::Solver::PseudoTransient::getSolverStatistics() const
{ return globalDataPtr->getSolverStatistics(); }
