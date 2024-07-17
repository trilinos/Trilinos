// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Observer_Vector.hpp"

void NOX::ObserverVector::runPreIterate(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPreIterate(solver);
}

void NOX::ObserverVector::runPostIterate(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPostIterate(solver);
}

void NOX::ObserverVector::runPreSolve(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPreSolve(solver);
}

void NOX::ObserverVector::runPostSolve(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPostSolve(solver);
}

void NOX::ObserverVector::runPreSolutionUpdate(const NOX::Abstract::Vector& update,
                                               const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPreSolutionUpdate(update,solver);
}

void NOX::ObserverVector::runPostSolutionUpdate(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPostSolutionUpdate(solver);
}

void NOX::ObserverVector::runPreLineSearch(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPreLineSearch(solver);
}

void NOX::ObserverVector::runPostLineSearch(const NOX::Solver::Generic& solver)
{
  for (auto& i : vec_)
    i->runPostLineSearch(solver);
}

void NOX::ObserverVector::pushBack(const Teuchos::RCP<NOX::Observer>& observer)
{
  vec_.push_back(observer);
}

void NOX::ObserverVector::popBack()
{
  vec_.pop_back();
}

void NOX::ObserverVector::clear()
{
  vec_.clear();
}
