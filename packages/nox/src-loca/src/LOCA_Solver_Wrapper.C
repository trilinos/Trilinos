// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#include "LOCA_Solver_Wrapper.H"	// class definition
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

NOX::StatusTest::StatusType 
LOCA::Solver::Wrapper::getStatus()
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

const NOX::Abstract::Group& 
LOCA::Solver::Wrapper::getPreviousSolutionGroup() const
{
  return *oldSolnGrpPtr;
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

void 
LOCA::Solver::Wrapper::resetWrapper()
{
  // Get current and old solution groups
  const NOX::Abstract::Group& soln = constSolverPtr->getSolutionGroup();
  const NOX::Abstract::Group& oldSoln = 
    constSolverPtr->getPreviousSolutionGroup();

  const LOCA::Extended::MultiAbstractGroup* eGrpPtr;
  const LOCA::Extended::MultiAbstractGroup* oldEGrpPtr;

  // Cast soln group to an extended group
  eGrpPtr = dynamic_cast<const LOCA::Extended::MultiAbstractGroup*>(&soln);

  if (eGrpPtr == NULL) {
    // soln group is not extended, so set points to original groups
    solnGrpPtr = Teuchos::rcp(&soln, false);
    oldSolnGrpPtr = Teuchos::rcp(&oldSoln, false);
  }

  else {
    // soln group is extended so get underlying groups
    oldEGrpPtr = 
      dynamic_cast<const LOCA::Extended::MultiAbstractGroup*>(&oldSoln);
    solnGrpPtr = eGrpPtr->getUnderlyingGroup();
    oldSolnGrpPtr = oldEGrpPtr->getUnderlyingGroup();
  }

  return;
}
