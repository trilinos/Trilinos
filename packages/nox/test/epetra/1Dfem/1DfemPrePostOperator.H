// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEST_1DFEM_PREPOSTOPERATOR_H
#define TEST_1DFEM_PREPOSTOPERATOR_H

#include "NOX_Common.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Utils.H"

class UserPrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  //! Ctor.
  UserPrePostOperator(const NOX::Utils& u);

  //! Destructor.
  ~UserPrePostOperator();

  void runPreIterate(const NOX::Solver::Generic& solver);

  void runPostIterate(const NOX::Solver::Generic& solver);

  void runPreSolve(const NOX::Solver::Generic& solver);

  void runPostSolve(const NOX::Solver::Generic& solver);

  void runPreLineSearch(const NOX::Solver::Generic& solver);

  void runPostLineSearch(const NOX::Solver::Generic& solver);

  int getNumRunPreIterate() const { return numRunPreIterate; };

  int getNumRunPostIterate() const { return numRunPostIterate; };

  int getNumRunPreSolve() const { return numRunPreSolve; };

  int getNumRunPostSolve() const { return numRunPostSolve; };

  int getNumRunPreLineSearch() const { return numRunPreLineSearch; };

  int getNumRunPostLineSearch() const { return numRunPostLineSearch; };

protected:

  NOX::Utils utils;

  int numRunPreIterate;
  int numRunPostIterate;
  int numRunPreSolve;
  int numRunPostSolve;
  int numRunPreLineSearch;
  int numRunPostLineSearch;

};
#endif
