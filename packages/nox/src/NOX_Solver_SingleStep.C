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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Solver_SingleStep.H"    // class definition
#include "NOX_GlobalData.H"    // class definition
#include "NOX_Abstract_Group.H"    // class definition
#include "NOX_Abstract_Group.H"    // class definition
#include "NOX_Solver_SolverUtils.H"
#include "Teuchos_ParameterList.hpp"  // class data element

NOX::Solver::SingleStep::
SingleStep(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
       const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(xGrp),                         // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  paramsPtr(p),
  ignoreLinearSolverFailures(false),
  updateJacobian(true)
{
  Teuchos::ParameterList validParams;
  validParams.set("Ignore Linear Solver Failures",false);
  validParams.set("Update Jacobian",true);
  p->sublist("Single Step Solver").validateParametersAndSetDefaults(validParams,0);
  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils();
  prePostOperator.reset(utilsPtr,p->sublist("Solver Options"));
  ignoreLinearSolverFailures = p->sublist("Single Step Solver").get<bool>("Ignore Linear Solver Failures");
  updateJacobian = p->sublist("Single Step Solver").get<bool>("Update Jacobian");
  if (not updateJacobian)
    frozenJacobianPtr = solnPtr->clone(DeepCopy); // take ownership of Jacobian
  init();
}

// Protected
void NOX::Solver::SingleStep::init()
{
  // Initialize
  nIter = 0;
  status = NOX::StatusTest::Unconverged;

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

}

void NOX::Solver::SingleStep::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  init();
}

void NOX::Solver::SingleStep::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>&)
{
  solnPtr->setX(initialGuess);
  init();
}

NOX::Solver::SingleStep::~SingleStep()
{

}

NOX::StatusTest::StatusType NOX::Solver::SingleStep::getStatus()
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

  if (updateJacobian) {
    if (!check(solnPtr->computeJacobian(), "compute Jacobian"))
      return false;
  }

  // Pick Jacobian matrix to use
  Teuchos::RCP<NOX::Abstract::Group> jacobian = updateJacobian ? solnPtr : frozenJacobianPtr;

  // Reuse memory in group instead of new allocation for dir
  NOX::Abstract::Vector& dir = const_cast<NOX::Abstract::Vector&>(solnPtr->getNewton());
  const auto status = jacobian->applyJacobianInverse(paramsPtr->sublist("Linear Solver"),
                                                     solnPtr->getF(),
                                                     dir);

  if (!ignoreLinearSolverFailures) {
    if (!check(status,"solve Newton system"))
      return false;
  }

  solnPtr->computeX(*oldSolnPtr, dir, -1.0);

  return true;
}

NOX::StatusTest::StatusType NOX::Solver::SingleStep::step()
{
  prePostOperator.runPreIterate(*this);

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;

  if (try_step())
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Failed;

  prePostOperator.runPostIterate(*this);

  printUpdate();

  return status;
}

NOX::StatusTest::StatusType NOX::Solver::SingleStep::solve()
{
  prePostOperator.runPreSolve(*this);

  step();

  prePostOperator.runPostSolve(*this);

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

// protected
void NOX::Solver::SingleStep::printUpdate()
{
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- The \"Nonlinear\" Solver Step -- \n";
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }
}
