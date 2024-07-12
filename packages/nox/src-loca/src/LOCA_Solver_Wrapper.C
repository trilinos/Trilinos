// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Solver_Wrapper.H"    // class definition
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_Extended_MultiAbstractGroup.H"

LOCA::Solver::Wrapper::
Wrapper(const Teuchos::RCP<NOX::Solver::Generic>& solver) :
  solverPtr(solver),
  constSolverPtr(solver)
{
  resetWrapper();
}

LOCA::Solver::Wrapper::
Wrapper(const Teuchos::RCP<const NOX::Solver::Generic>& solver) :
  solverPtr(Teuchos::null),
  constSolverPtr(solver)
{
  resetWrapper();
}

LOCA::Solver::Wrapper::~Wrapper()
{
}

void
LOCA::Solver::Wrapper::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& tests)
{
  solverPtr->reset(initialGuess, tests);
  resetWrapper();
}

void
LOCA::Solver::Wrapper::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solverPtr->reset(initialGuess);
  resetWrapper();
}

void
LOCA::Solver::Wrapper::
reset()
{
  resetWrapper();
}

NOX::StatusTest::StatusType
LOCA::Solver::Wrapper::getStatus() const
{
  return solverPtr->getStatus();
}

NOX::StatusTest::StatusType
LOCA::Solver::Wrapper::step()
{
  NOX::StatusTest::StatusType status = solverPtr->step();
  resetWrapper();
  return status;
}

NOX::StatusTest::StatusType
LOCA::Solver::Wrapper::solve()
{
  NOX::StatusTest::StatusType status = solverPtr->solve();
  resetWrapper();
  return status;
}

const NOX::Abstract::Group&
LOCA::Solver::Wrapper::getSolutionGroup() const
{
  return *solnGrpPtr;
}

Teuchos::RCP< const NOX::Abstract::Group >
LOCA::Solver::Wrapper::getSolutionGroupPtr() const
{
  return solnGrpPtr;
}

const NOX::Abstract::Group&
LOCA::Solver::Wrapper::getPreviousSolutionGroup() const
{
  return *oldSolnGrpPtr;
}

Teuchos::RCP< const NOX::Abstract::Group >
LOCA::Solver::Wrapper::getPreviousSolutionGroupPtr() const
{
  return oldSolnGrpPtr;
}

int
LOCA::Solver::Wrapper::getNumIterations() const
{
  return constSolverPtr->getNumIterations();
}

const Teuchos::ParameterList&
LOCA::Solver::Wrapper::getList() const
{
  return constSolverPtr->getList();
}

Teuchos::RCP< const Teuchos::ParameterList >
LOCA::Solver::Wrapper::getListPtr() const
{
  return constSolverPtr->getListPtr();
}

Teuchos::RCP<const NOX::SolverStats>
LOCA::Solver::Wrapper::getSolverStatistics() const
{
  return constSolverPtr->getSolverStatistics();
}

void
LOCA::Solver::Wrapper::resetWrapper()
{
  // Get current and old solution groups
  const NOX::Abstract::Group& soln = constSolverPtr->getSolutionGroup();
  const NOX::Abstract::Group& oldSoln =
    constSolverPtr->getPreviousSolutionGroup();

  const LOCA::Extended::MultiAbstractGroup* eGrpPtr;

  // Cast soln group to an extended group
  eGrpPtr = dynamic_cast<const LOCA::Extended::MultiAbstractGroup*>(&soln);

  if (eGrpPtr == NULL) {
    // soln group is not extended, so set points to original groups
    solnGrpPtr = Teuchos::rcp(&soln, false);
    oldSolnGrpPtr = Teuchos::rcp(&oldSoln, false);
  } else {
    // soln group is extended so get underlying groups
    const LOCA::Extended::MultiAbstractGroup & oldEGrp =
      dynamic_cast<const LOCA::Extended::MultiAbstractGroup&>(oldSoln);
    solnGrpPtr = eGrpPtr->getUnderlyingGroup();
    oldSolnGrpPtr = oldEGrp.getUnderlyingGroup();
  }

  return;
}
