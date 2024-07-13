// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "1DfemPrePostOperator.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Solver_Generic.H"


UserPrePostOperator::UserPrePostOperator(const NOX::Utils& u) :
  numRunPreIterate(0),
  numRunPostIterate(0),
  numRunPreSolve(0),
  numRunPostSolve(0),
  numRunPreLineSearch(0),
  numRunPostLineSearch(0)
{
  utils = u;
}

UserPrePostOperator::~UserPrePostOperator()
{

}

void UserPrePostOperator::
runPreIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPreIterate;
  utils.out(NOX::Utils::Details) <<
    "1Dfem's runPreIterate() routine called!" << std::endl;
}

void UserPrePostOperator::
runPostIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPostIterate;
  utils.out(NOX::Utils::Details)
    << "1Dfem's runPostIterate() routine called!" << std::endl;
}

void UserPrePostOperator::
runPreSolve(const NOX::Solver::Generic& solver)
{
  ++numRunPreSolve;
  utils.out(NOX::Utils::Details)
    << "1Dfem's runPreSolve() routine called!" << std::endl;
}

void UserPrePostOperator::
runPostSolve(const NOX::Solver::Generic& solver)
{
  ++numRunPostSolve;
  utils.out(NOX::Utils::Details)
    << "1Dfem's runPostSolve() routine called!" << std::endl;
}

void UserPrePostOperator::
runPreLineSearch(const NOX::Solver::Generic& solver)
{
  ++numRunPreLineSearch;
  utils.out(NOX::Utils::Details)
    << "1Dfem's runPreLineSearch() routine called!" << std::endl;
}

void UserPrePostOperator::
runPostLineSearch(const NOX::Solver::Generic& solver)
{
  ++numRunPostLineSearch;
  utils.out(NOX::Utils::Details)
    << "1Dfem's runPostLineSearch() routine called!" << std::endl;
}
