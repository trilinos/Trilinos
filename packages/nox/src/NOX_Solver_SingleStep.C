// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Solver_SingleStep.H"    // class definition
#include "NOX_GlobalData.H"    // class definition
#include "NOX_Abstract_Group.H"    // class definition
#include "NOX_Abstract_Group.H"    // class definition
#include "NOX_Solver_SolverUtils.H"
#include "Teuchos_ParameterList.hpp"  // class data element
#include "NOX_Observer.hpp"
#include "NOX_SolverStats.hpp"

NOX::Solver::SingleStep::
SingleStep(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
       const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(xGrp),                         // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  paramsPtr(p),
  ignoreLinearSolverFailures(false),
  updateJacobian(true),
  printNorms(false),
  computeRelativeNorm(false),
  normF_0(0.0)
{
  Teuchos::ParameterList validParams;
  validParams.set("Ignore Linear Solver Failures",false,"Return the step as converged, ignoring the returned linear solver status.");
  validParams.set("Update Jacobian",true,"Recompute the Jacobian at each Newton iteration.");
  validParams.set("Print Norms",false,"Print the norms at each iteration.");
  validParams.set("Compute Relative Norm",false, "Computes relative norm, print only if \"Print Norms\" is enabled.");
  validParams.sublist("Linear Solver"); // Allows for arbitrary/user defined linear solve parameters to be set
  p->sublist("Single Step Solver").validateParametersAndSetDefaults(validParams,0);
  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils();
  observer = NOX::Solver::parseObserver(p->sublist("Solver Options"));
  ignoreLinearSolverFailures = p->sublist("Single Step Solver").get<bool>("Ignore Linear Solver Failures");
  updateJacobian = p->sublist("Single Step Solver").get<bool>("Update Jacobian");
  if (!updateJacobian)
    frozenJacobianPtr = solnPtr->clone(DeepCopy); // take ownership of Jacobian
  printNorms =  p->sublist("Single Step Solver").get<bool>("Print Norms");
  computeRelativeNorm =  p->sublist("Single Step Solver").get<bool>("Compute Relative Norm");
  init();
}

// Protected
void NOX::Solver::SingleStep::init()
{
  // Initialize
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

}

void NOX::Solver::SingleStep::
reset()
{
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();
}

void NOX::Solver::SingleStep::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();
}

void NOX::Solver::SingleStep::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>&)
{
  solnPtr->setX(initialGuess);
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  globalDataPtr->getNonConstSolverStatistics()->reset();
}

NOX::Solver::SingleStep::~SingleStep()
{

}

NOX::StatusTest::StatusType NOX::Solver::SingleStep::getStatus() const
{
  return status;
}

bool NOX::Solver::SingleStep::check(NOX::Abstract::Group::ReturnType ret,
    const std::string& task)
{
  if (ret != NOX::Abstract::Group::Ok) {
    if (utilsPtr->isPrintType(Utils::Error))
      utilsPtr->out() << "NOX::Solver::SingleStep - Unable to " << task << std::endl;
    return false;
  }
  return true;
}

bool NOX::Solver::SingleStep::try_step()
{
  if (!check(solnPtr->computeF(), "compute F"))
    return false;

  if (computeRelativeNorm) {
    normF_0 = solnPtr->getF().norm();
  }

  if (updateJacobian) {
    if (!check(solnPtr->computeJacobian(), "compute Jacobian"))
      return false;
  }

  // Pick Jacobian matrix to use
  Teuchos::RCP<NOX::Abstract::Group> jacobian = updateJacobian ? solnPtr : frozenJacobianPtr;

  // Reuse memory in group instead of new allocation for dir
  NOX::Abstract::Vector& dir = const_cast<NOX::Abstract::Vector&>(solnPtr->getNewton());
  const auto ls_status = jacobian->applyJacobianInverse(paramsPtr->sublist("Single Step Solver").sublist("Linear Solver"),
                                                        solnPtr->getF(),
                                                        dir);

  jacobian->logLastLinearSolveStats(*globalDataPtr->getNonConstSolverStatistics());

  if (!ignoreLinearSolverFailures) {
    if (!check(ls_status,"solve Newton system"))
      return false;
  }

  observer->runPreSolutionUpdate(dir,*this);
  solnPtr->computeX(*oldSolnPtr, dir, -1.0);
  observer->runPostSolutionUpdate(*this);

  return true;
}

NOX::StatusTest::StatusType NOX::Solver::SingleStep::step()
{
  observer->runPreIterate(*this);

  // SingleStep solver means step() is always a new solve
  globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearSolves();

  // Update iteration count.
  nIter ++;
  globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearIterations();

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;

  if (try_step())
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Failed;

  observer->runPostIterate(*this);

  printUpdate();

  return status;
}

NOX::StatusTest::StatusType NOX::Solver::SingleStep::solve()
{
  observer->runPreSolve(*this);

  this->reset();

  step();

  observer->runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group&
NOX::Solver::SingleStep::getSolutionGroup() const
{
  return *solnPtr;
}

Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::SingleStep::getSolutionGroupPtr() const
{
  return solnPtr;
}

const NOX::Abstract::Group&
NOX::Solver::SingleStep::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::SingleStep::getPreviousSolutionGroupPtr() const
{
  return oldSolnPtr;
}

int NOX::Solver::SingleStep::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList&
NOX::Solver::SingleStep::getList() const
{
  return *paramsPtr;
}

Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::SingleStep::getListPtr() const
{
   return paramsPtr;
}

Teuchos::RCP<const NOX::SolverStats>
NOX::Solver::SingleStep::getSolverStatistics() const
{
  return globalDataPtr->getSolverStatistics();
}

// protected
void NOX::Solver::SingleStep::printUpdate()
{
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- The \"Nonlinear\" Solver Step -- \n";
    if (printNorms) {
      if (!solnPtr->isF())
        solnPtr->computeF();
      double normF = solnPtr->getF().norm();
      double normDx = solnPtr->getNewtonPtr()->norm();
      utilsPtr->out() << "||F||=" << normF << ", ||dx||=" << normDx;
      if (computeRelativeNorm) {
        utilsPtr->out() << ", ||F|| / ||F_0||=" << normF/normF_0;
      }
    }
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }
}
