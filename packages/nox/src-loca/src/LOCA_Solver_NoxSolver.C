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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Solver_NoxSolver.H"	  // class definition

//LOCA includes
#include "LOCA_Abstract_Group.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Utils.H"

// NOX includes
#include "NOX_Solver_Manager.H"
#include "NOX_Common.H"

using namespace LOCA;
using namespace LOCA::Solver;

NoxSolver::NoxSolver(Abstract::Group& xGrp, 
		     NOX::StatusTest::Generic& t, 
		     const NOX::Parameter::List& p) :
  Generic(p),
  initialParams(p),
  statusTest(t),
  initialGroupPtr(0),
  solverGroup(xGrp),
  noxSolverPtr(0),
  label("NOX Nonlinear Solver")
{
  initialGroupPtr = xGrp.clone(NOX::DeepCopy);
  noxSolverPtr = new NOX::Solver::Manager(xGrp, t, p.sublist("Solver"));
}

NoxSolver::~NoxSolver() 
{
  delete initialGroupPtr;
  delete noxSolverPtr;
}

bool NoxSolver::reset(ResetType t)
{
  if (t == InitialGuess) {
    solverGroup = *initialGroupPtr;
  }
  return noxSolverPtr->reset(solverGroup, statusTest, initialParams);
}

bool NoxSolver::setParams(ParameterVector& p)
{
  const LOCA::Abstract::Group& grp = dynamic_cast<const LOCA::Abstract::Group&>(noxSolverPtr->getSolutionGroup());

  return (const_cast<LOCA::Abstract::Group&>(grp)).setParams(p);
}

NOX::StatusTest::StatusType NoxSolver::getStatus()
{
  return noxSolverPtr->getStatus();
}

NOX::StatusTest::StatusType NoxSolver::iterate()
{
  // Call computeF - nox expects it to be calculated
  (const_cast<LOCA::Abstract::Group&>(getSolutionGroup())).computeF();

  return noxSolverPtr->iterate();
}

NOX::StatusTest::StatusType NoxSolver::solve()
{
  // Call computeF - nox expects it to be calculated
  (const_cast<LOCA::Abstract::Group&>(getSolutionGroup())).computeF();

  return noxSolverPtr->solve();
}

const Abstract::Group& NoxSolver::getSolutionGroup() const
{
  const LOCA::Abstract::Group& locaGroup = dynamic_cast<const LOCA::Abstract::Group&>(noxSolverPtr->getSolutionGroup());
  return locaGroup;
}

const Abstract::Group& NoxSolver::getPreviousSolutionGroup() const
{
  const LOCA::Abstract::Group& locaGroup = dynamic_cast<const LOCA::Abstract::Group&>(noxSolverPtr->getPreviousSolutionGroup());
  return locaGroup;
}

int NoxSolver::getNumIterations() const
{
  return noxSolverPtr->getNumIterations();
}

NOX::Parameter::List& NoxSolver::getParameterList()
{
  // NOX duplicates the parameter list passed into the Solver::Manager ctor
  // so we have to copy that into the loca parameter list
  const NOX::Parameter::List& noxParamList = noxSolverPtr->getParameterList();
  params.sublist("Solver") = noxParamList;
  return params;
}

const NOX::Parameter::List& NoxSolver::getParameterList() const
{
  // NOX duplicates the parameter list passed into the Solver::Manager ctor
  // so we have to copy that into the loca parameter list
  const NOX::Parameter::List& noxParamList = noxSolverPtr->getParameterList();
  params.sublist("Solver") = noxParamList;
  return params;
}

const string NoxSolver::getLabel() const
{
  return label;
}
