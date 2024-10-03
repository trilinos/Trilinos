// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_THYRA_OBSERVER_PRINT_TEST_HPP
#define NOX_THYRA_OBSERVER_PRINT_TEST_HPP

#include "NOX_Common.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Utils.H"

class ObserverPrintTest : public NOX::Abstract::PrePostOperator {

public:

  //! Ctor.
  ObserverPrintTest(const NOX::Utils& u);

  //! Destructor.
  ~ObserverPrintTest();

  void runPreIterate(const NOX::Solver::Generic& solver);

  void runPostIterate(const NOX::Solver::Generic& solver);

  void runPreSolve(const NOX::Solver::Generic& solver);

  void runPostSolve(const NOX::Solver::Generic& solver);
  
  int getNumPreIterateCalls() const {return numPreIterateCalls;}
  int getNumPostIterateCalls() const {return numPostIterateCalls;}
  int getNumPreSolveCalls() const {return numPreSolveCalls;}
  int getNumPostSolveCalls() const {return numPostSolveCalls;}

protected:

  NOX::Utils utils;

  int numPreIterateCalls;
  int numPostIterateCalls;
  int numPreSolveCalls;
  int numPostSolveCalls;

};
#endif
