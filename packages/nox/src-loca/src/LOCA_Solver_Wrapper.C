// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Solver_Wrapper.H"	// class definition
#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_Extended_AbstractGroup.H"

LOCA::Solver::Wrapper::
Wrapper(NOX::Solver::Generic& solver) :
  solverPtr(&solver)
{
  resetWrapper();
}

LOCA::Solver::Wrapper::
Wrapper(const NOX::Solver::Generic& solver) :
  solverPtr(&(const_cast<NOX::Solver::Generic&>(solver)))
{
  resetWrapper();
}

LOCA::Solver::Wrapper::~Wrapper()
{
}

bool
LOCA::Solver::Wrapper::
reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& params)
{
  bool res = solverPtr->reset(grp, tests, params);
  resetWrapper();
  return res;
}

bool 
LOCA::Solver::Wrapper::
reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests) 
{
  bool res =  solverPtr->reset(grp, tests);
  resetWrapper();
  return res;
}

NOX::StatusTest::StatusType 
LOCA::Solver::Wrapper::getStatus()
{
  return solverPtr->getStatus();
}

NOX::StatusTest::StatusType 
LOCA::Solver::Wrapper::iterate()
{
  NOX::StatusTest::StatusType status = solverPtr->iterate();
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
  return solverPtr->getNumIterations();
}

const Teuchos::ParameterList& 
LOCA::Solver::Wrapper::getList() const
{
  return solverPtr->getList();
}

void 
LOCA::Solver::Wrapper::resetWrapper()
{
  // Get current and old solution groups
  const NOX::Abstract::Group& soln = solverPtr->getSolutionGroup();
  const NOX::Abstract::Group& oldSoln = solverPtr->getPreviousSolutionGroup();

  const LOCA::Extended::AbstractGroup* eGrpPtr;
  const LOCA::Extended::AbstractGroup* oldEGrpPtr;

  // Cast soln group to an extended group
  eGrpPtr = dynamic_cast<const LOCA::Extended::AbstractGroup*>(&soln);

  if (eGrpPtr == NULL) {
    // soln group is not extended, so set points to original groups
    solnGrpPtr = Teuchos::rcp(&soln, false);
    oldSolnGrpPtr = Teuchos::rcp(&oldSoln, false);
  }

  else {
    // soln group is extended so get underlying groups
    oldEGrpPtr = dynamic_cast<const LOCA::Extended::AbstractGroup*>(&oldSoln);
    solnGrpPtr = Teuchos::rcp(&(eGrpPtr->getUnderlyingGroup()), false);
    oldSolnGrpPtr = Teuchos::rcp(&(oldEGrpPtr->getUnderlyingGroup()), false);
  }

  return;
}
