// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
